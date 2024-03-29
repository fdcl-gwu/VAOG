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


T=10;

[t_GD, J_GD, errR_GD] = gdSO3(A,J,R0,T);
figure(1);plot(t_GD ,J_GD, 'k');
set(gca,'xscale','log','yscale','log');
hold on;
figure(2);plot(t_GD, errR_GD, 'k');
set(gca,'xscale','log','yscale','log');
hold on;


p=2;

[t_AGD, J_AGD, errR_AGD] = agdSO3(A,p,J,R0,T);
figure(1);plot(t_AGD, J_AGD, 'b:');
figure(2);plot(t_AGD, errR_AGD, 'b:');

hk = 0.1;
flag_update_h = true;
[t_LGVI, J_LGVI, errR_LGVI] = vigdSO3(A,p,J,R0,hk,T,flag_update_h);

figure(1);plot(t_LGVI, J_LGVI, 'r:');
figure(2);plot(t_LGVI, errR_LGVI, 'r:');

p=4;

[t_AGD, J_AGD, errR_AGD] = agdSO3(A,p,J,R0,T);
figure(1);plot(t_AGD, J_AGD, 'b');
figure(2);plot(t_AGD, errR_AGD, 'b');

hk = 0.1;
flag_update_h = true;
[t_LGVI, J_LGVI, errR_LGVI] = vigdSO3(A,p,J,R0,hk,T,flag_update_h);

figure(1);plot(t_LGVI, J_LGVI, 'r');
figure(2);plot(t_LGVI, errR_LGVI, 'r');

p=6;

[t_AGD, J_AGD, errR_AGD] = agdSO3(A,p,J,R0,T);
figure(1);plot(t_AGD, J_AGD, 'b', 'linewidth', 2);
figure(2);plot(t_AGD, errR_AGD, 'b','linewidth', 2);

hk = 0.1;
flag_update_h = true;
[t_LGVI, J_LGVI, errR_LGVI] = vigdSO3(A,p,J,R0,hk,T,flag_update_h);

figure(1);plot(t_LGVI, J_LGVI, 'r', 'linewidth', 2);
figure(2);plot(t_LGVI, errR_LGVI, 'r', 'linewidth', 2);


grid on;

save('comp_0');

figure(1);
xlabel('$t$','interpreter','latex');
ylabel('$f-f^*$','interpreter','latex');

figure(2);
xlabel('$t$','interpreter','latex');
ylabel('$\|I-R^TR\|$','interpreter','latex');
