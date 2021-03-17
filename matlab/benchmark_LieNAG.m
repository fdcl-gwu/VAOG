clc;
clear;
close all;

%% setup
A = [    0.7572    0.5678    0.5308
    0.7537    0.0759    0.7792
    0.3804    0.0540    0.9340];
r =     [0.7915
    -0.3101
    0.5267];
[U S V]=psvd(A);
R0 = U*expm(0.9*pi*hat(r))*V';


RIC=R0;
xiIC=zeros(3);

gamma=1;            % this needs to be tuned
h=0.01;              % this needs to be tuned
TotalSteps=10000;

%% Lie-NAG-SC and Lie-NAG-C
%   implementation: 2nd-order splitting (explicit and conformally symplectic)
%   it preserves Lie group to machine precision; optimization error also has this lower bound

for iter=1:2    % 1: NAG-SC, 2: NAG-C
    tic

    R=RIC;
    xi=xiIC;
    loss=zeros(1,TotalSteps+1);
    deviation=zeros(1,TotalSteps+1);
    i=0;    loss(i+1)=norm(A-R,'fro')^2/2;  deviation(i+1)=norm(R'*R-eye(3),'fro');
    RR=R;    xxi=xi;

    hh=h/2;
    temp=expm(hh*xxi);
    RR=RR*temp;
    
    for i=1:TotalSteps

        hh=h;
        force=-(A'*RR-RR'*A);
        if (iter==1)    % NAG-SC
            xxi=exp(-gamma*hh)*xxi+(1-exp(-gamma*hh))/gamma*force;
        end
        if (iter==2)    % NAG-C  gamma=3/((i-1/2)*h)
            t0=(i-1)*h; t1=i*h;
            xxi=t0^3/t1^3*xxi+(t1^4-t0^4)/t1^3/4*force;
        end

        hh=h;
        temp=expm(hh*xxi);
        RR=RR*temp;

        R=RR;
        xi=xxi;

        loss(i+1)=norm(A-R,'fro')^2/2;
        deviation(i+1)=norm(R'*R-eye(3),'fro');
    end

    switch iter
        case 1
            h_NAGSC=h;
            loss_NAGSC=loss;
            deviation_NAGSC=deviation;
        case 2
            h_NAGC=h;
            loss_NAGC=loss;
            deviation_NAGC=deviation;
    end

    toc
end

%%
%truth
[U,S,V]=psvd(A);
Rstar=U*V';
exact=norm(A-Rstar,'fro')^2/2;

figure
% subplot(2,1,1);
plot([0:TotalSteps],loss_NAGSC-exact, [0:TotalSteps],loss_NAGC-exact, 'LineWidth',2);
set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('iteration steps');
title('distance from optimum value');
% 
% subplot(2,1,2);
% semilogy([0:TotalSteps],deviation_NAGSC, [0:TotalSteps],deviation_NAGC, 'LineWidth',2);
% xlabel('iteration steps');
% title('deviation from Lie group');
% axis([0 TotalSteps 1E-16 1E-9]);

legend(['NAG-SC h=',num2str(h_NAGSC)], ['NAG-C h=',num2str(h_NAGC)], 'Location','Best');
