clear;
close all;

R0=expmso3(10*pi/180*[1 0.5 -0.5]');
x0 = [2, 1, -0.5]';

p=3;
J=eye(3);
h=0.0015;
N=5000;
Cx=1.5;

flag_update_h = false;
flag_display_progress = false;
flag_stop = 2;

[t, JJ, errR, t_elapsed, R, x] = vigdSE3(p,J,R0,x0, h,N,flag_update_h,flag_display_progress,flag_stop,Cx);
disp(JJ(end));
figure();
plot(1:N,JJ);
set(gca,'xscale','log','yscale','log');
title(['C_x =' num2str(Cx)]);
drawnow;



R=R(:,:,end);
x=x(:,end);


save comp_6

