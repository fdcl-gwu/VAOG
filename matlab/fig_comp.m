clear all;
close all;

load gdSO3_0
plot(t_GD ,J_GD, 'r', t_AGD, J_AGD,'b', t_AGD_Pi, J_AGD_Pi,'g');
set(gca,'xscale','log','yscale','log');
hold on;

load vigdSO3_0
plot(t,JJ,'k:');
