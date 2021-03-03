clear all;
close all;

p = rand(3,1);
J = rand(3,3);
J = J + J';
K=1/2*trace(inv(J))*eye(3)-inv(J);
R = expmso3(rand(3,1));
mu = rand(3,1);
F = expmso3(rand(3,1));
P = R*F*hat(mu);

% 1/2*mu'*inv(J)*mu
% 1/2*trace(inv(J)*mu*mu')
% -1/2*trace(K*hat(mu)^2)
% 
% return;

Q=R'*P;
eta = rand(3,1);
h=rand(1);

trace(-hat(eta)*R'*P)/2 - h/2*trace(K*-hat(eta)*Q*Q') - h/2*trace(K*Q*Q'*hat(eta))
trace(-hat(eta)*Q)/2 + h*trace(hat(eta)*Q*Q'*K) 
eta'*vee( (Q - Q')/2 - h*Q*Q'*K + h*K*Q*Q')
eta'*vee( (Q - Q')/2 + h*hat(F*mu)^2*K - h*K*hat(F*mu)^2)
eta'* ( (trace(F)*eye(3)-F')/2 + h*hat(K*F*mu)*F)*mu


delP = rand(3,3);
1/2*trace(R'*delP) - h*trace(K*R'*delP*P'*R)
1/2*trace(delP'*( R - 2*h*R*K*R'*P))



Rkp = expmso3(rand(3,1));
delp = rand(3,1);
delP = Rkp*hat(delp)
F = R'*Rkp;

1/2*trace(delP'*( -Rkp + R - 2*h*R*K*R'*P))
1/2*trace(hat(delp)'*( -eye(3) + F' - 2*h*F'*K*R'*P))


