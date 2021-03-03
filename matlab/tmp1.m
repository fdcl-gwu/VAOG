clear all;
close all;


J = rand(3,3);
J = J + J';
K = 1/2*trace(inv(J))*eye(3)-inv(J);
Jd = 1/2*trace(J)*eye(3) - J;

F = expmso3(0.05*rand(3,1))

%mu = vee(Jd*F-F'*Jd)
mu = vee(F*Jd-Jd*F')

2*(trace(K)*eye(3)-K)*mu
vee(F-F')