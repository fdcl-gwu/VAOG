function [JJ, errR, t_elapsed] = benchmark_LieNAG_new(A,p,R0,hk,N,flag_iter)

RIC=R0;
xiIC=zeros(3);

C=1;
CC=C*p^2;  % this needs to be tuned
gamma=1;            % this needs to be tuned
h=hk;              % this needs to be tuned
TotalSteps=N;

%% Lie-NAG-SC, Lie-NAG-C, and Lie-Bregman
%   implementation: 2nd-order splitting (explicit and conformally symplectic)
%   it preserves Lie group to machine precision; optimization error also has this lower bound

% 1: NAG-SC, 2: NAG-C, 3: Bregman

iter=flag_iter;


R=RIC;
xi=xiIC;
loss=zeros(1,TotalSteps);
deviation=zeros(1,TotalSteps);

RR=R;    xxi=xi;
RR_save=zeros(3,3,N);
RR_save(:,:,1)=RR;

tic
hh=h/2;
temp=expm(hh*xxi);
RR=RR*temp;

for i=1:TotalSteps-1
    
    hh=h;       % note: i'm using the same h for all 3 methods, but Bregman actually requires smaller h, especially if t is large.
    % early stopping could be one option to sell Bregman better
    force=-(A'*RR-RR'*A);
    switch iter
        case 1    % NAG-SC
            xxi=exp(-gamma*hh)*xxi+(1-exp(-gamma*hh))/gamma*force;
        case 2    % NAG-C
            t0=(i-1)*h; t1=i*h;
            xxi=t0^3/t1^3*xxi+(t1^4-t0^4)/t1^3/4*force;
        case 3    % more general Bregman (p \neq 2)
            t0=(i-1)*h; t1=i*h;
            xxi=t0^(p+1)/t1^(p+1)*xxi+(t1^(2*p)-t0^(2*p))/t1^(1+p)/(2*p)*CC*force;
    end
    
    hh=h;
    
    temp=expmso3(hh*vee(xxi));
    RR=RR*temp;
    
    R=RR;
    xi=xxi;
    RR_save(:,:,i+1)=RR;
end

t_elapsed=toc;

for i=1:TotalSteps
    loss(i)=norm(A-RR_save(:,:,i),'fro')^2/2;
    deviation(i)=norm(RR_save(:,:,i)'*RR_save(:,:,i)-eye(3),'fro');
end

%%
%truth
[U,S,V]=svd(A);
Rstar=U*diag([1 1 det(U*V)])*V';
exact=norm(A-Rstar,'fro')^2/2;

JJ = loss-exact;
errR = deviation;
