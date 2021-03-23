function tmp_comp_grad

A=rand(3,3);



R=expmso3(rand(3,1));

eps = 1e-8;
M = grad(R,A);
M_num = grad_num(R,A,eps);

disp(norm(M_num-M));

end


function J = obj(R,A)
tmp=R-A;
J = tmp(:)'*tmp(:)/2;
end

function M = grad(R,A)
% this actually computes the negative gradient
M=vee(R'*A-A'*R);
end

function M = grad_num(R,A,eps)
% this actually computes the negative gradient
eps = 1e-8;

II=eye(3);
M=zeros(3,1);
for i=1:3
    M(i)=(obj(R*expmso3(eps*II(:,i)),A)-obj(R,A))/eps;
end
M=-M;

end



