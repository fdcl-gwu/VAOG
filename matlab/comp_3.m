% effects of varying h_0
clear all;
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

N = 10000;
flag_update_h = false;
flag_display_progress = false;
flag_stop = 2;
p = 3;

h0=0.005;
[t, JJ, errR, t_elapsed, N_grad_comp] = vigdSO3(A,p,J,R0,h0,N,flag_update_h,flag_display_progress,flag_stop);
disp(t_elapsed);

figure(1);plot(1:length(t) ,JJ);
set(gca,'xscale','log','yscale','log');

hold on;

for flag_type=1:3
    [JJ, errR, t_elapsed] = benchmark_LieNAG_new(A,p,R0,h0,N,flag_type);
    figure(1);plot(1:length(JJ),JJ,'--');
    disp(t_elapsed);
end


save comp_3


% figure(1);plot(t ,JJ, 'k');
% set(gca,'xscale','log','yscale','log');
% hold on;
% figure(2);plot(t, errR, 'k');
% set(gca,'xscale','log','yscale','log');
% hold on;
% 
