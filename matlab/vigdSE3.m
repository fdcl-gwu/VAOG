function [t, JJ, errR, t_elapsed, R, x] = vigdSE3(p,J,R0,x0,hk,TN,flag_update_h,flag_display_progress,flag_stop,Cx)
% flag_stop =1 : stop of t>TN,
% flag_stop =2 : stop of k>TN,
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
P = [p_W_landmarks'; ones(1,N)];
p_new = K_normal*K*[eye(3), zeros(3,1)]*P;
keypoints_n = p_new(1:2,:)./p_new(3,:);
keypoints_n = [keypoints_n(2,:); keypoints_n(1,:)];
K_n = K_normal*K;

keypoints=keypoints_n;
K=K_n;

JJ_opt = obj(eye(3),zeros(3,1), K,keypoints,p_W_landmarks);
Jd=trace(J)/2*eye(3)-J;%nonstandard inertia matrix

options = optimoptions(@lsqnonlin,'StepTolerance', 1e-4, 'Display', 'off');

C=1;

%%
if flag_stop == 2
    R=zeros(3,3,TN);
    x=zeros(3,TN);
    v=zeros(3,TN);
    grad_x=zeros(3,TN);
    grad_R=zeros(3,TN);
    t=zeros(TN,1);
    Pi=zeros(3,TN);
end

R(:,:,1)=R0;
Pi(:,1)=zeros(3,1);
x(:,1)=x0;
v(:,1)=zeros(3,1);

t(1)=1;
JJ(1) = obj(R(:,:,1),x(:,1), K,keypoints,p_W_landmarks);
[grad_R(:,1), grad_x(:,1)] = grad(R(:,:,1),x(:,1), K,keypoints,p_W_landmarks);

tic;
if flag_update_h
%     [~, E(1)] = err_FLdm(hk, t(1), Pi(:,1), R(:,:,1), 0, JJ(1), grad_R(:,1), p, C, Jd, K,keypoints,p_W_landmarks);
end
k = 1;


while 1
    
    if flag_update_h
%         [hk, res] = lsqnonlin(@(hk) err_FLdm(hk, t(k), Pi(:,k), R(:,:,k), E(k), JJ(k), grad_R(:,k), p, C, Jd, K,keypoints,p_W_landmarks), hk, 0, [], options);
    end
    t(k+1) = t(k) + hk;
    tkkp= (t(k)+t(k+1))/2;

    
    [phi_kkp, phi_dot_kkp] = compute_phi(tkkp,p);
    [theta_k, ~] = compute_theta(t(k),p,C);
    [theta_kp, theta_dot_kp] = compute_theta(t(k+1),p,C);

  
    g = (Pi(:,k) - hk*theta_k/2*grad_R(:,k))*hk/phi_kkp;    
    if norm(g) > 1
        disp('Warning: g is too large');
    end
    
    Fk = expmso3( g*asin(norm(g))/norm(g) );    
    R(:,:,k+1) = R(:,:,k)*Fk;
    
    delxk = (v(:,k) - hk*theta_k/2*Cx*grad_x(:,k))*hk/phi_kkp;
    x(:,k+1) = x(:,k)+delxk;    

    [grad_R(:,k+1), grad_x(:,k+1)] = grad(R(:,:,k+1),x(:,k+1), K,keypoints,p_W_landmarks);       

    v(:,k+1) = v(:,k) - hk*theta_k/2*Cx*grad_x(:,k) - hk*theta_kp/2*Cx*grad_x(:,k+1);    
    Pi(:,k+1) = Fk'*Pi(:,k) - hk*theta_k/2*Fk'*grad_R(:,k) - hk*theta_kp/2*grad_R(:,k+1);   
    
    JJ(k+1) = obj(R(:,:,k+1),x(:,k+1), K,keypoints,p_W_landmarks);
%     if flag_update_h
%         JJ(k+1) = obj(R(:,:,k+1),K,keypoints,p_W_landmarks);
% 
%         E(k+1) = (-phi_dot_kkp/2/hk+phi_kkp/hk^2)*trace((eye(3)-Fk)*Jd) ...
%             + hk*theta_dot_kp/2*JJ(k+1) + theta_k/2*JJ(k) +theta_kp/2*JJ(k+1);
%         
%     end
    
    k=k+1;
    
    switch flag_stop
        case 1
            if t(k) > TN
                break;
            end
            if flag_display_progress && rem(k,500)==0
                disp([t(k)/TN]);
            end
            
        case 2
            if k >= TN
                break;
            end
            if flag_display_progress && rem(k,500)==0
                disp([k/TN]);
            end
    end
    
end
t_elapsed=toc;

%%

for k=1:length(t)
%     JJ(k) = obj(R(:,:,k),x(:,k), K,keypoints,p_W_landmarks)-JJ_opt;
    errR(k) = norm(R(:,:,k)'*R(:,:,k)-eye(3));
end

end

function [errE FLdm]= err_FLdm(hk, tk, Pik, Rk, Ek, JJk, grad_Rk, p, C, Jd, K,keypoints,p_W_landmarks)
tkp = tk + hk;
tkkp= (tk+tkp)/2;

[theta_k, theta_dot_k] = compute_theta(tk,p,C);
[theta_kp, ~] = compute_theta(tkp,p,C);
[phi_kkp, phi_dot_kkp] = compute_phi(tkkp,p);

g = (Pik - hk/2*theta_k*grad_Rk)*hk/phi_kkp;
Fk = expmso3( g*asin(norm(g))/norm(g) );
Rkp = Rk*Fk;

JJkp = obj(Rkp,K,keypoints,p_W_landmarks);

FLdm = (phi_dot_kkp/2/hk + phi_kkp/hk^2)*trace((eye(3)-Fk)*Jd) ...
    -hk*theta_dot_k/2*JJk + theta_k/2*JJk + theta_kp/2*JJkp;

errE = abs(Ek-FLdm);
end

function J = obj(R, x, K,keypoints,p_W_landmarks)
N = length(keypoints);
p = [keypoints([2,1],:); ones(1,N)];
P = [p_W_landmarks'; ones(1,N)];

M = K*[R, x];
p_proj = M*P;

J = norm(p(1:2,:)./p(3,:)-p_proj(1:2,:)./p_proj(3,:))^2;

end

function [grad_R, grad_x] = grad(R, x, K,keypoints,p_W_landmarks)
eps = 1e-8;
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


N = length(keypoints);
p = [keypoints([2,1],:); ones(1,N)];
P = [p_W_landmarks'; ones(1,N)];

M = K*[R, x];
p_proj = M*P;
uv = p(1:2,:)./p(3,:);
uv_proj = p_proj(1:2,:)./p_proj(3,:);

J = norm(uv-uv_proj)^2/N/N;

grad_x=zeros(3,1);

for i=1:N
    tmp = 2*(uv_proj(:,i)-uv(:,i))';
    grad_x(1) = grad_x(1) + tmp*[K(1,1)/(x(3) + P(1:3,i)'*R(3,:)'); 0];
    grad_x(2) = grad_x(2) + tmp*[0; K(1,1)/(x(3) + P(1:3,i)'*R(3,:)')];
    grad_x(3) = grad_x(3) + tmp*[-(K(1,1)*(x(1) + P(1:3,i)'*R(1,:)'))/(x(3) + P(1:3,i)'*R(3,:)')^2;...
        -(K(1,1)*(x(2) + P(1:3,i)'*R(2,:)'))/(x(3) + P(1:3,i)'*R(3,:)')^2];
end
% grad_x=grad_x/N/N;


end

function [phi, phi_dot] = compute_phi(t,p)
phi = t^(p+1)/p;
phi_dot = (p+1)/p*t^p;
end

function [theta, theta_dot] = compute_theta(t,p,C)
theta = C*p*t^(2*p-1);
theta_dot = C*p*(2*p-1)*t^(2*p-2);
end




