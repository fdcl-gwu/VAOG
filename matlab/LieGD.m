function [JJ, errR, t_elapsed] = LieGD(A,R0,h,N)
% flag_stop =1 : stop of t>TN,
% flag_stop =2 : stop of k>TN,

[U S V]=psvd(A);
R_opt = U*V';
JJ_opt = obj(R_opt,A);

%%
R(:,:,1)=R0;


tic;
for k=1:N-1  
    Mk = grad(R(:,:,k),A);    
    R(:,:,k+1) = R(:,:,k)*expmso3(h*Mk);
    
end
t_elapsed=toc;

%%

for k=1:N
    JJ(k) = obj(R(:,:,k),A)-JJ_opt;
    errR(k) = norm(R(:,:,k)'*R(:,:,k)-eye(3));
end

end

function J = obj(R,A)
tmp=R-A;
J = tmp(:)'*tmp(:)/2;
end

function M = grad(R,A)
M=vee(R'*A-A'*R);
end


