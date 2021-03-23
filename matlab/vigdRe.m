function [t, JJ, t_elapsed, x] = vigdRe(p,x0,hk,TN,flag_update_h,flag_display_progress,flag_stop)
% flag_stop =1 : stop of t>TN,
% flag_stop =2 : stop of k>TN,

n=length(x0);

A=rand(n,n);
A=A+A';
[V D]=eig(A);
A= V*abs(D)*V';
disp(eig(A));

b=rand(n,1);

x_opt = -A\b;
JJ_opt = obj(x_opt,A,b);


options = optimoptions(@lsqnonlin,'StepTolerance', 1e-12,'OptimalityTolerance', 1e-6, 'Display', 'off');
% options = optimoptions(@lsqnonlin,'StepTolerance', 1e-8);

C=1;
%%
if flag_stop == 2
    x=zeros(n,TN);
    grad_x=zeros(n,TN);
    t=zeros(TN,1);
    v=zeros(n,TN);
end

x(:,1)=x0;
v(:,1)=zeros(n,1);
t(1)=1;
JJ(1) = obj(x(:,1),A,b);
grad_x(:,1) = grad(x(:,1),A,b);

tic;
if flag_update_h
    [~, E(1)] = err_FLdm(hk, t(1), v(:,1), x(:,1), 0, JJ(1), grad_x(:,1), p, C, A,b);    
end
k = 1;


while 1
  
    if flag_update_h
        [hk, res] = lsqnonlin(@(hk) err_FLdm(hk, t(k), v(:,k), x(:,k), E(k), JJ(k), grad_x(:,1), p, C, A,b), 1, 0, [], options);
    end
%     disp(hk);
    t(k+1) = t(k) + hk;
    tkkp= (t(k)+t(k+1))/2;

    [phi_kkp, phi_dot_kkp] = compute_phi(tkkp,p);
    [theta_k, ~] = compute_theta(t(k),p,C);
    [theta_kp, theta_dot_kp] = compute_theta(t(k+1),p,C);
    
    delxk = (v(:,k)-hk*theta_k/2*grad_x(:,k))*hk/phi_kkp;
    x(:,k+1) = x(:,k)+delxk;    

    grad_x(:,k+1) = grad(x(:,k+1),A,b);
    v(:,k+1) = v(:,k) - hk*theta_k/2*grad_x(:,k) - hk*theta_kp/2*grad_x(:,k+1);
    
    JJ(k+1) = obj(x(:,k+1),A,b);
    if flag_update_h
        
        E(k+1) = (-phi_dot_kkp/2/hk+phi_kkp/hk^2)*norm(delxk)^2/2 ...
            + hk*theta_dot_kp/2*JJ(k+1) + theta_k/2*JJ(k) +theta_kp/2*JJ(k+1);
    end
    
    k=k+1;
    
    switch flag_stop
        case 1
            if t(k) > TN
                break;
            end
            if flag_display_progress && rem(k,500)==0
                disp([t(k)/TN, hk]);
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


JJ = JJ-JJ_opt;
% for k=1:length(t)
%     JJ(k) = obj(x(:,k),A,b)-JJ_opt;
%     errR(k) = norm(R(:,:,k)'*R(:,:,k)-eye(3));
% end



end

function [errE, FLdm]= err_FLdm(hk, tk, vk, xk, Ek, JJk, grad_k, p, C, A, b)
tkp = tk + hk;
tkkp= (tk+tkp)/2;

[theta_k, theta_dot_k] = compute_theta(tk,p,C);
[theta_kp, ~] = compute_theta(tkp,p,C);
[phi_kkp, phi_dot_kkp] = compute_phi(tkkp,p);

delxk = (vk-hk*theta_k/2*grad_k)*hk/phi_kkp;
xkp = xk+delxk;
JJkp = obj(xkp,A,b);

FLdm = (phi_dot_kkp/2/hk + phi_kkp/hk^2)*norm(delxk)^2/2 ...
    -hk*theta_dot_k/2*JJk + theta_k/2*JJk + theta_kp/2*JJkp;

errE = abs(Ek-FLdm);
end

function [phi, phi_dot] = compute_phi(t,p)
phi = t^(p+1)/p;
phi_dot = (p+1)/p*t^p;
end

function [theta, theta_dot] = compute_theta(t,p,C)
theta = C*p*t^(2*p-1);
theta_dot = C*p*(2*p-1)*t^(2*p-2);
end

function J = obj(x,A,b)
J=1/2*(A*x+b)'*(A*x+b);
end

% function M = grad(R,A)
% % this actually computes the negative gradient
% M=vee(R'*A-A'*R);
% end

function M = grad(x,A,b)
M = A'*A*x+A'*b;
end



