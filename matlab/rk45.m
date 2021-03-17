function [t, JJ, errR, t_elapsed, sol] = rk45(A,p,J,R0,T)

[U S V]=psvd(A);
R_opt = U*V';
JJ_opt = obj(R_opt,A);

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

C=1;

%%
W0 = zeros(3,1);

X0= [reshape(R0,9,1); W0];

tic;
sol = ode45(@(t,X) eom_EL_Breg(t,X,A,p,J,C), [1, T+1], X0, options);
t_elapsed=toc;

t=sol.x;X=sol.y;

for k=1:length(t)
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


