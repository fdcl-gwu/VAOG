function [t_AGD, J_AGD, errR_AGD]=agdSO3(A,p,J,R0,T)
global A


% A=rand(3,3);%diag([5, 3, 1])*expm(pi*hat([1 1 1])/sqrt(3));
% A = [    0.7572    0.5678    0.5308
%     0.7537    0.0759    0.7792
%     0.3804    0.0540    0.9340];

[U S V]=psvd(A);
R_opt = U*V';
J_opt = obj(R_opt);

% p=4;
C=1;
% J=eye(3);

%% AGD
W0=zeros(3,1);
X0=[reshape(R0,9,1); W0];

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
sol=ode45(@(t,X) eom_AGD(t,X,p,J,C), [1 T], X0, options);
t=sol.x;X=sol.y';disp(sol.stats);

N=length(t);
W=X(:,10:12)';

for k=1:N
    R_AGD(:,:,k) = reshape(X(k,1:9),3,3);
    J_AGD(k) = obj(R_AGD(:,:,k))-J_opt;
    Pi_AGD(:,k) = t(k)^(p+1)/p*J*W(:,k);
    errR_AGD(k) = norm(R_AGD(:,:,k)'*R_AGD(:,:,k)-eye(3));
end
t_AGD = t;

% %% AGD_Pi
% 
% 
% X0=[reshape(R0,9,1); J*W0];
% 
% sol=ode45(@eom_AGD_Pi, tspan, X0, options);
% t=sol.x;X=sol.y';disp(sol.stats);
% 
% N=length(t);
% Pi=X(:,10:12)';
% 
% for k=1:N
%     R_AGD_Pi(:,:,k) = reshape(X(k,1:9),3,3);
%     J_AGD_Pi(k) = obj(R_AGD_Pi(:,:,k))-J_opt;
% end
% t_AGD_Pi = t;
% 
% 
% %% Post Process
% 
% plot(t_GD ,J_GD, 'r', t_AGD, J_AGD,'b', t_AGD_Pi, J_AGD_Pi,'g');
% set(gca,'xscale','log','yscale','log');

% filename='gdSO3_0';
% save(filename);
% evalin('base','clear all');
% evalin('base',['load ' filename]);

end

function J = obj(R)
global A

tmp=R-A;
J = tmp(:)'*tmp(:)/2;

end

function M = grad(R)
global A

M=vee(R'*A-A'*R);%+randn(3,1)*1e-3;

end


function X_dot = eom_GD(t,X)
R=reshape(X,3,3);
M = grad(R);

R_dot = R*hat(M);
X_dot = reshape(R_dot,9,1);
end

function X_dot = eom_AGD(t,X, p, J,C)
R=reshape(X(1:9),3,3);
W=X(10:12);

M = grad(R);
R_dot = R*hat(W);
W_dot = J\( -(p+1)/t*J*W - hat(W)*J*W + C*p^2*t^(p-2)*M);

X_dot = [reshape(R_dot,9,1); W_dot];
end

function X_dot = eom_AGD_Pi(t,X)
global p J C

R=reshape(X(1:9),3,3);

Pi=X(10:12);

M = grad(R);
W = J\(p/t^(p+1)*Pi);

R_dot = R*hat(W);
Pi_dot = -hat(W)*Pi + C*p*t^(2*p-1)*M;

X_dot = [reshape(R_dot,9,1); Pi_dot];
end
