clear all;
close all;

load lgviSO3 h J R W

K = 1/2*trace(inv(J))*eye(3) - inv(J);

R1 = R(:,:,1);
mu2 = J*W(:,2);

mu_kp = mu2;

delF = 1;
F = eye(3);
while delF > 1e-15
    F_new = findF( 2*h*(trace(K)*eye(3) - F'*K*F)*mu_kp, 2*eye(3));   
    delF = norm(F_new-F);
    F = F_new;
end

R2 = R1*F-R(:,:,2)

mu_k = (1/2*(trace(F)*eye(3)-F') - h*hat(K*F*mu_kp)*F)*mu_kp;
mu_k - J*W(:,1)