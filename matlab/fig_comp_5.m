% checking the numerical gradient evaluation
clear;
close all;

% A = [    0.7572    0.5678    0.5308
%     0.7537    0.0759    0.7792
%     0.3804    0.0540    0.9340];
% r =     [0.7915
%     -0.3101
%     0.5267];

% A=rand(3,3);
% r=rand(3,1);


A =[
    0.9398    0.6393    0.5439
    0.6456    0.5447    0.7210
    0.4795    0.6473    0.5225];

r =[
    0.9937
    0.2187
    0.1058
];

[U S V]=psvd(A);
R0 = U*expm(pi*hat(r))*V';
J=eye(3);

N = 10000;


%% LGVI
flag_update_h = false;
flag_display_progress = false;
flag_stop = 2;


p = 2;        h = 0.1;
p = 3;        h = 0.025;

[t, JJ, errR, t_elapsed] = vigdSO3(A,p,J,R0,h,N,flag_update_h,flag_display_progress,flag_stop);
disp(t_elapsed);

figure(1);plot(1:N,JJ);
hold on;
set(gca,'xscale','log','yscale','log');


p = 4;        h = 0.0055;
[t, JJ, errR, t_elapsed] = vigdSO3(A,p,J,R0,h,N,flag_update_h,flag_display_progress,flag_stop);
disp(t_elapsed);

figure(1);plot(1:N,JJ);


figure(1)
ylabel('${f}(R)-{f}(R^*)$','interpreter','latex');
xlabel('$k$','interpreter','latex');
set(gca,'FontSize',16);

