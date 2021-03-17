function [t, JJ, errR, t_elapsed] = rk2(A,p,J,R0,h,N)

[U S V]=psvd(A);
R_opt = U*V';
JJ_opt = obj(R_opt,A);

Jd=trace(J)/2*eye(3)-J;%nonstandard inertia matrix

options = optimoptions(@lsqnonlin,'StepTolerance', 1e-4, 'Display', 'off');

N_grad_comp = 0;
C=1;

t=1:h:(N-1)*h+1;
%%

W0 = zeros(3,1);

X=zeros(12,N);
X(:,1)= [reshape(R0,9,1); W0];
% options = odeset('RelTol',1e-8,'AbsTol',1e-8);
tic;

for k=1:N-1
  F1 = eom_EL_Breg(t(k),X(:,k),A,p,J,C);
  F2 = eom_EL_Breg(t(k)+h/2,X(:,k)+h/2*F1,A,p,J,C);
  F3 = eom_EL_Breg(t(k)+h/2,X(:,k)+h/2*F2,A,p,J,C);
  F4 = eom_EL_Breg(t(k)+h,X(:,k)+h*F3,A,p,J,C);
  
  X(:,k+1) = X(:,k) + (h/6)*(F1 + 2*F2 + 2*F3 + F4);

end
t_elapsed=toc;


for k=1:N
    R(:,:,k)=reshape(X(1:9,k),3,3);
    JJ(k) = obj(R(:,:,k),A)-JJ_opt;
    errR(k) = norm(R(:,:,k)'*R(:,:,k)-eye(3));
end

end

function X_dot = eom_EL_Breg(t,X,A,p,J,C)
R=reshape(X(1:9),3,3);
W=X(10:12);

M = grad(R,A);
R_dot = R*hat(W);
W_dot = J\( -(p+1)/t*J*W - hat(W)*J*W + C*p^2*t^(p-2)*M);

X_dot = [reshape(R_dot,9,1); W_dot];
end


function J = obj(R,A)

tmp=R-A;
J = tmp(:)'*tmp(:)/2;

end

function M = grad(R,A)
global N_grad_comp

N_grad_comp = N_grad_comp+1;
M=vee(R'*A-A'*R);
end


