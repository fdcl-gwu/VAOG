function [t, JJ, errR, t_elapsed] = vigdSO3(A,p,J,R0,hk,TN,flag_update_h,flag_display_progress,flag_stop)
% flag_stop =1 : stop of t>TN,
% flag_stop =2 : stop of k>TN,

[U S V]=psvd(A);
R_opt = U*V';
JJ_opt = obj(R_opt,A);
Jd=trace(J)/2*eye(3)-J;%nonstandard inertia matrix

options = optimoptions(@lsqnonlin,'StepTolerance', 1e-4, 'Display', 'off');

C=1;
%%
if flag_stop == 2
    R=zeros(3,3,TN);
    M=zeros(3,TN);
    t=zeros(TN,1);
    Pi=zeros(3,TN);
end

R(:,:,1)=R0;
Pi(:,1)=zeros(3,1);
t(1)=1;
JJ(1) = obj(R(:,:,1),A);
M(:,1) = grad(R(:,:,1),A);

tic;
[~, E(1)] = err_FLdm(hk, t(1), Pi(:,1), R(:,:,1), 0, JJ(1), M(:,1), p, C, Jd, A);
k = 1;

while 1
    
    %     if flag_update_h
    %         [hk, res] = lsqnonlin(@(hk) err_FLdm(hk, t(k), Pi(:,k), R(:,:,k), E(k), JJ(k), M(:,k), p, C, Jd, A), hk, 0, [], options);
    %     end
    t(k+1) = t(k) + hk;
    tkkp= (t(k)+t(k+1))/2;
    g = (Pi(:,k) + hk/2*C*p*t(k)^(2*p-1)*M(:,k))*hk*p/tkkp^(p+1);
    Fk = expmso3( g*asin(norm(g))/norm(g) );
    
    R(:,:,k+1) = R(:,:,k)*Fk;
    M(:,k+1) = grad(R(:,:,k+1),A);
    Pi(:,k+1) = Fk'*Pi(:,k) + hk/2*C*p*t(k)^(2*p-1)*Fk'*M(:,k) + hk/2*C*p*t(k+1)^(2*p-1)*M(:,k+1);
    
    
    %     if flag_update_h
    %         JJ(k+1) = obj(R(:,:,k+1),A);
    %         E(k+1) = (-(p+1)*tkkp^p/2/hk/p + tkkp^(p+1)/hk^2/p)*trace((eye(3)-Fk)*Jd)...
    %             +1/2*C*p*t(k)^(2*p-1)*JJ(k) + 1/2*C*p*t(k+1)^(2*p-1)*JJ(k+1) ...
    %             +hk/2*C*p*(2*p-1)*t(k+1)^(2*p-2)*JJ(k+1);
    %     end
    
    k=k+1;
    
    switch flag_stop
        case 1
            if t(k) > TN
                break;
            end
            %             if flag_display_progress && rem(k,500)==0
            %                 disp([t(k)/TN]);
            %             end
            
        case 2
            if k >= TN
                break;
            end
            % %             if flag_display_progress && rem(k,500)==0
            % %                 disp([k/TN]);
            %             end
    end
    
end
t_elapsed=toc;

%%

for k=1:length(t)
    JJ(k) = obj(R(:,:,k),A)-JJ_opt;
    errR(k) = norm(R(:,:,k)'*R(:,:,k)-eye(3));
end

end

function [errE FLdm]= err_FLdm(hk, tk, Pik, Rk, Ek, Jk, Mk, p, C, Jd, A)
tkp = tk + hk;
tkkp= (tk+tkp)/2;
g = (Pik + hk/2*C*p*tk^(2*p-1)*Mk)*hk*p/tkkp^(p+1);

Fk = expmso3( g*asin(norm(g))/norm(g) );
Rkp = Rk*Fk;
Jkp = obj(Rkp,A);

FLdm = ((p+1)*tkkp^p/2/hk/p + tkkp^(p+1)/hk^2/p)*trace((eye(3)-Fk)*Jd)...
    +1/2*C*p*tk^(2*p-1)*Jk + 1/2*C*p*tkp^(2*p-1)*Jkp ...
    -hk/2*C*p*(2*p-1)*tk^(2*p-2)*Jk;

errE = abs(Ek-FLdm);
end


function J = obj(R,A)
tmp=R-A;
J = tmp(:)'*tmp(:)/2;
end

function M = grad(R,A)
M=vee(R'*A-A'*R);
end


