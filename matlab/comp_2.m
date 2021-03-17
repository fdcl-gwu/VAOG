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

T = 30
flag_update_h = true;
flag_display_progress = true;
flag_stop =1;
hk = 0.1;
p = 4;

h0_vec=[0.001 0.05 0.01 0.1 0.4];

for i=1:length(h0_vec)
    h0 = h0_vec(i)
    [t, JJ, errR, t_elapsed, N_grad_comp] = vigdSO3(A,p,J,R0,h0,T,flag_update_h,flag_display_progress,flag_stop);
    t_vec{i} = t;
    J_vec{i} = JJ;
    errR_vec{i} = errR;
    t_elapsed_vec{i} = t_elapsed;
    N_grad_comp_vec{i} = N_grad_comp;
    disp([hk log10(JJ(end)) log10(errR(end)) t_elapsed N_grad_comp(end)]);
    hold on;
end

save comp_2


% figure(1);plot(t ,JJ, 'k');
% set(gca,'xscale','log','yscale','log');
% hold on;
% figure(2);plot(t, errR, 'k');
% set(gca,'xscale','log','yscale','log');
% hold on;
% 
