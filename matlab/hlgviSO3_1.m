clear all;
close all;

load lgviSO3 t h J R W N E

K = 1/2*trace(inv(J))*eye(3) - inv(J);

Q(:,:,1) = R(:,:,1);
mu(:,1) = J*W(:,1);
eps = 1e-12;

for k=1:N-1   
    delF = 1;
    Fk = eye(3);
    while delF > eps
        gk = 2*h*(trace(K)*eye(3)-K)*Fk'*mu(:,k);
        Fk_new = findF(gk, 2*eye(3));
        delF = norm(Fk_new -Fk);
        Fk = Fk_new;
    end
    
    F(:,:,k)= Fk;
    mu(:,k+1) = F(:,:,k)'*mu(:,k);
    Q(:,:,k+1) = Q(:,:,k)*F(:,:,k);

    disp(norm( vee(Fk-Fk') - 2*h*(trace(K)*eye(3)-K)*mu(:,k+1) ));
end

for k=1:N
    Wh(:,k) = inv(J)*mu(:,k);
%     H(k) = 1/2*Wh(:,k)'*J*Wh(:,k);
    H(k) = -1/2*trace(K*hat(mu(:,k))^2);
    errR(k) = norm(R(:,:,k)-Q(:,:,k));
end

plot(t,W,'r',t,Wh,'b:')
figure;
plot(t,E,'r',t,H, 'b:')

