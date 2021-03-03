clear all;
close all;

load lgviSO3 t h J R W N

K = 1/2*trace(inv(J))*eye(3) - inv(J);
FF = @(f, mu) (eye(3) - hat(f) - f*f.') * inv(2*h*(trace(K)*eye(3) - K)) * 4*f / (1+f.'*f)^2 - mu;
caley_f2F = @(f) (eye(3)-hat(f))\(eye(3)+hat(f));
caley_F2f = @(F) vee(inv(F+eye(3))*(F-eye(3)));

Q(:,:,1) = R(:,:,1);
mu(:,1) = J*W(:,1);

options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12);

for k=1:N-1
    fk = fsolve(@(f) FF(f,mu(:,k)), zeros(3,1));
    F(:,:,k) = caley_f2F(fk);
    Q(:,:,k+1) = Q(:,:,k)*F(:,:,k);
    mu(:,k+1) = inv(2*h*(trace(K)*eye(3) - K))*vee(F(:,:,k)-F(:,:,k)');
end

for k=1:N
    Wh(:,k) = inv(J)*mu(:,k);
    H(k) = 1/2*Wh(:,k)'*J*Wh(:,k);
    errR(k) = norm(R(:,:,k)-Q(:,:,k));
end

% plot(t,W,'r',t,Wh,'b:')
% figure;
% plot(t,H)

for k=1:N-1
    err1(k) = norm( mu(:,k) - 1/2*vee( -F(:,:,k)*hat(mu(:,k+1))  +hat(mu(:,k+1))*F(:,:,k)') );
    err11(k) = norm( mu(:,k) - 1/2*(trace(F(:,:,k)*eye(3)-F(:,:,k)'))*mu(:,k+1) );
    err2(k) = norm( 2*h*(trace(K)*eye(3)-K)*mu(:,k+1) - vee(F(:,:,k)-F(:,:,k)') );
end
