% comparision with other discretization of Bregman Euler-Lagrange equation
clear;
close all;

A = [    0.7572    0.5678    0.5308
    0.7537    0.0759    0.7792
    0.3804    0.0540    0.9340];
r =     [0.7915
    -0.3101
    0.5267];
[U S V]=psvd(A);
R0 = U*expm(0.9*pi*hat(r))*V';
J=eye(3);

N = 9000;
p = 6;
h = 0.001;

%% LGVI
flag_update_h = false;
flag_display_progress = false;
flag_stop = 2;
[t, JJ, errR, t_elapsed] = vigdSO3(A,p,J,R0,h,N,flag_update_h,flag_display_progress,flag_stop);
disp(t_elapsed);

figure(1);plot(t,JJ);
hold on;
set(gca,'xscale','log','yscale','log');
figure(2);plot(t, errR);
set(gca,'xscale','log','yscale','log');
hold on;

%% Splitting
flag_type=3;
[JJ, errR, t_elapsed] = benchmark_LieNAG_new(A,p,R0,h,N,flag_type);
figure(1);plot(t,JJ,'--');
figure(2);plot(t, errR, '--');
disp(t_elapsed);

%% RK4
[t, JJ, errR, t_elapsed] = rk4(A,p,J,R0,h,N);
figure(1);plot(t ,JJ, '--');
figure(2);plot(t, errR, '--');
disp(t_elapsed);

%% RK45
[t, JJ, errR, t_elapsed, sol] = rk45(A,p,J,R0,h*N);
figure(1);plot(t ,JJ, '-.');
figure(2);plot(t, errR, '-.');

disp(t_elapsed);
save comp_3

%%
figure(1)
legend('$\mathrm{ELGVI}$', '$\mathrm{Splitting (SPLT)}$', '$\mathrm{RK4}$', '$\mathrm{RK45}$',...
    'Location','SouthWest','interpreter','latex')
ylabel('${f}(R)-{f}(R^*)$','interpreter','latex');
xlabel('$t$','interpreter','latex');
set(gca,'FontSize',16);

figure(2)
legend('$\mathrm{ELGVI}$', '$\mathrm{Splitting\; (SPLT)}$', '$\mathrm{RK4}$', '$\mathrm{RK45}$',...
    'Location','SouthEast','interpreter','latex')

xlabel('$t$','interpreter','latex');
ylabel('$\|I-R^TR\|$','interpreter','latex');
set(gca,'FontSize',16);

figure(1);
print('comp_3a','-depsc');
figure(2);
print('comp_3b','-depsc');

cpu_time = [    0.0774    0.0648    0.0919    0.1074    0.1185    0.1454    0.0868    0.0788    0.0711    0.0854
    0.0231    0.0235    0.0274    0.0262    0.0245    0.0379    0.0230    0.0264    0.0219    0.0243
    0.3346    0.3899    0.3970    0.4047    0.3916    0.3796    0.3972    0.4212    0.3748    0.3562
    1.0050    1.1003    1.2049    1.1693    1.1611    1.1086    1.1300    1.2166    1.2657    1.1144];


disp(mean(cpu_time'))
  
    
!cp comp_3?.eps ../figs
