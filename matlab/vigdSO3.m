function vigdSO3
%% second order attempt

global A p J C Jd
close all

% A=rand(3,3);%diag([5, 3, 1])*expm(pi*hat([1 1 1])/sqrt(3));
% A = eye(3);
% A = [ 0.6536    0.1765    0.5259
%     0.4579    0.8653    0.9180
%     0.8375    0.9539    0.3723];

load gdSO3_0 A

[U S V]=psvd(A);
R_opt = U*V';
J_opt = obj(R_opt);

p=3;
C=1;
J=eye(3);
Jd=trace(J)/2*eye(3)-J;%nonstandard inertia matrix


%%
T = 30;
hk = 0.6;

R(:,:,1)=eye(3);
Pi(:,1)=zeros(3,1);
t(1)=1;
JJ(1) = obj(R(:,:,1));
M(:,1) = grad(R(:,:,1));

[~, E(1)] = err_FLdm(hk, t(1), Pi(:,1), R(:,:,1), 0, JJ(1), M(:,1));

options = optimoptions(@lsqnonlin,'StepTolerance', 1e-8, 'Display', 'off');

k = 1;
while t(k) < T
    if rem(k,100) == 0
        disp([t(k)/T]);
    end
    
    [hk, res] = lsqnonlin(@(hk) err_FLdm(hk, t(k), Pi(:,k), R(:,:,k), E(k), JJ(k), M(:,k)), 0.01, 0, [], options);
    
    t(k+1) = t(k) + hk;
    tkkp= (t(k)+t(k+1))/2;
    g = (Pi(:,k) + hk/2*C*p*t(k)^(2*p-1)*M(:,k))*hk*p/tkkp^(p+1);
    Fk = expmso3( g*asin(norm(g))/norm(g) );
    
    R(:,:,k+1) = R(:,:,k)*Fk;
    JJ(k+1) = obj(R(:,:,k+1));
    M(:,k+1) = grad(R(:,:,k+1));
    Pi(:,k+1) = Fk'*Pi(:,k) + hk/2*C*p*t(k)^(2*p-1)*Fk'*M(:,k) + hk/2*C*p*t(k+1)^(2*p-1)*M(:,k+1);
    E(k+1) = (-(p+1)*tkkp^p/2/hk/p + tkkp^(p+1)/hk^2/p)*trace((eye(3)-Fk)*Jd)...
        +1/2*C*p*t(k)^(2*p-1)*JJ(k) + 1/2*C*p*t(k+1)^(2*p-1)*JJ(k+1) ...
        +hk/2*C*p*(2*p-1)*t(k+1)^(2*p-2)*JJ(k+1);
    
    k=k+1;
end

%%

JJ=JJ-J_opt;
plot(t,JJ);
set(gca,'xscale','log','yscale','log');


filename='vigdSO3_0';
save(filename);
evalin('base','clear all');
evalin('base',['load ' filename]);

end

function [errE FLdm]= err_FLdm(hk, tk, Pik, Rk, Ek, Jk, Mk)
global p C Jd

tkp = tk + hk;
tkkp= (tk+tkp)/2;
g = (Pik + hk/2*C*p*tk^(2*p-1)*Mk)*hk*p/tkkp^(p+1);

Fk = expmso3( g*asin(norm(g))/norm(g) );
Rkp = Rk*Fk;
Jkp = obj(Rkp);

FLdm = ((p+1)*tkkp^p/2/hk/p + tkkp^(p+1)/hk^2/p)*trace((eye(3)-Fk)*Jd)...
    +1/2*C*p*tk^(2*p-1)*Jk + 1/2*C*p*tkp^(2*p-1)*Jkp ...
    -hk/2*C*p*(2*p-1)*tk^(2*p-2)*Jk;

errE = abs(Ek-FLdm);
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


