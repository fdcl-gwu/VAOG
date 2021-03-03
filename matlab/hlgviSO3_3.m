clear all;
close all;

load lgviSO3 t h J R W N E m g rho e3 Jd

Q=R;

R(:,:,1) = Q(:,:,1);
Pi(:,1) = J*W(:,1);
eps = 1e-10;
M(:,1)=m*g*cross(rho,R(:,:,1)'*e3);

% Pi2F = @(Pi) expm(h*hat(J\Pi));
% errPi = @(Pi,Y) (trace(Pi2F(Pi))*eye(3)-Pi2F(Pi)')*Pi - Y;

options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12);

tic;
for k=1:N-1
%     if rem(k,10) == 0
%         disp(k);
%     end
    
    F(:,:,k) = findF(h*Pi(:,k), J);
    R(:,:,k+1) = R(:,:,k)*F(:,:,k);    
    M(:,k+1)=m*g*cross(rho,R(:,:,k+1)'*e3);
    
    Yk = (trace(F(:,:,k)*Jd)*eye(3) - F(:,:,k)*Jd)'*inv(J)*Pi(:,k)+h*M(:,k+1);
    
%     Pikp = fsolve(@(Pi) errPi(Pi,Yk), Pi(:,k), options);
%   fixed point iteration is faster than fsolve
    delPi = 1;
    Pikp = Pi(:,:,k);
    while delPi > eps
        Fkp = findF(h*Pikp, J);
        Pikp_new = ((trace(Fkp*Jd)*eye(3)-Fkp*Jd)'*inv(J))\(Yk);
        delPi = norm(Pikp_new - Pikp);
        Pikp = Pikp_new;
    end
    Pi(:,:,k+1)=Pikp;
end
toc;

for k=1:N
    Wh(:,k) = inv(J)*Pi(:,k);
    H(k) = 1/2*Wh(:,k)'*J*Wh(:,k)-m*g*rho'*R(:,:,k)'*e3;
    errR(k) = norm(R(:,:,k)-Q(:,:,k));
end

figure;
subplot(3,1,1);
plot(t,W,'r',t,Wh,'b:')
subplot(3,1,2);
plot(t,E,'r',t,H,'b:')
subplot(3,1,3);
plot(t,errR)
