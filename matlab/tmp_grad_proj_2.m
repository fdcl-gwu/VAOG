function tmp_grad_proj_2

load data/p_W_landmarks.txt;
load data/K.txt;
load data/keypoints.txt;
keypoints = keypoints';


img.w = 1271;
img.h = 376;
K_normal = [2/img.w 0 -1;
    0 2/img.h -1
    0 0 1];

N = size(keypoints,2);
p = [keypoints(2,:); keypoints(1,:); ones(1,N)];
P = [p_W_landmarks'; ones(1,N)];

p_new = K_normal*K*[eye(3), zeros(3,1)]*P;
keypoints_n = p_new(1:2,:)./p_new(3,:);
keypoints_n = [keypoints_n(2,:); keypoints_n(1,:)];
K_n = K_normal*K;


R=eye(3);
x=[0 0 0]';

J = obj(R, x, K_n,keypoints_n,p_W_landmarks);



[num_grad_R, num_grad_x] = grad_num(R, x, K_n,keypoints_n,p_W_landmarks);

M = K_n*[R, x];
p_proj = M*P;
uv = p(1:2,:)./p(3,:);
uv_proj = p_proj(1:2,:)./p_proj(3,:);

K = K_n;
grad_x=zeros(3,1);
grad_R=zeros(3,1);
for i=1:N
    tmp = 2*(uv_proj(:,i)-uv(:,i))';
    grad_x(1) = grad_x(1) + tmp*[K(1,1)/(x(3) + P(1:3,i)'*R(3,:)'); 0];
    grad_x(2) = grad_x(2) + tmp*[0; K(1,1)/(x(3) + P(1:3,i)'*R(3,:)')];
    grad_x(3) = grad_x(3) + tmp*[-(K(1,1)*(x(1) + P(1:3,i)'*R(1,:)')); ...
        -(K(1,1)*(x(2) + P(1:3,i)'*R(2,:)'))]/(x(3) + P(1:3,i)'*R(3,:)')^2;
        
    grad_R(1) = grad_R(1) + tmp* ...
        [(K(1,1)*(P(2,i)*R(1,3)*x(3) - P(3,i)*R(1,2)*x(3) - P(2,i)*R(3,3)*x(1) + P(3,i)*R(3,2)*x(1)))/x(3)^2;
        (K(1,1)*(P(2,i)*R(2,3)*x(3) - P(3,i)*R(2,2)*x(3) - P(2,i)*R(3,3)*x(2) + P(3,i)*R(3,2)*x(2)))/x(3)^2];
    
    grad_R(2) = grad_R(2) + tmp* ...
        [-(K(1,1)*(P(1,i)*R(1,3)*x(3) - P(3,i)*R(1,1)*x(3) - P(1,i)*R(3,3)*x(1) + P(3,i)*R(3,1)*x(1)))/x(3)^2;
        -(K(1,1)*(P(1,i)*R(2,3)*x(3) - P(3,i)*R(2,1)*x(3) - P(1,i)*R(3,3)*x(2) + P(3,i)*R(3,1)*x(2)))/x(3)^2];

    grad_R(3) = grad_R(3) + tmp* ...
        [(K(1,1)*(P(1,i)*R(1,2)*x(3) - P(2,i)*R(1,1)*x(3) - P(1,i)*R(3,2)*x(1) + P(2,i)*R(3,1)*x(1)))/x(3)^2;
        (K(1,1)*(P(1,i)*R(2,2)*x(3) - P(2,i)*R(2,1)*x(3) - P(1,i)*R(3,2)*x(2) + P(2,i)*R(3,1)*x(2)))/x(3)^2];


end
grad_x=grad_x/N/N;
grad_R=grad_R/N/N;

disp([num_grad_x, grad_x]);
disp([num_grad_R, grad_R]);

end


function J = obj(R, x, K,keypoints,p_W_landmarks)
N = size(keypoints,2);
p = [keypoints(2,:); keypoints(1,:); ones(1,N)];
P = [p_W_landmarks'; ones(1,N)];

M = K*[R, x];
p_proj = M*P;

J = norm(p(1:2,:)./p(3,:)-p_proj(1:2,:)./p_proj(3,:));

end

function [grad_R, grad_x] = grad_num(R, x, K,keypoints,p_W_landmarks)
eps = 1e-6;
II=eye(3);

grad_R=zeros(3,1);
grad_x=zeros(3,1);

for i=1:3
    grad_R(i)=(obj(R*expmso3(eps*II(:,i)), x, K,keypoints,p_W_landmarks)...
        -obj(R*expmso3(-eps*II(:,i)), x, K,keypoints,p_W_landmarks))/eps/2;
    grad_x(i)=(obj(R,x+eps*II(:,i), K,keypoints,p_W_landmarks)...
        -obj(R,x-eps*II(:,i), K,keypoints,p_W_landmarks))/eps/2;
end

% for i=1:3
%     grad_R(i)=(obj(R*expmso3(eps*II(:,i)), x, K,keypoints,p_W_landmarks)...
%         -obj(R, x, K,keypoints,p_W_landmarks))/eps;
%     grad_x(i)=(obj(R,x+eps*II(:,i), K,keypoints,p_W_landmarks)...
%         -obj(R, x, K,keypoints,p_W_landmarks))/eps;
% end

end
