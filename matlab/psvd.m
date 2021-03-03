function [U S V]=psvd(F)

[U S V]=svd(F);
S=S*diag([1 1 det(U*V)]);
U=U*diag([1 1 det(U)]);
V=V*diag([1 1 det(V)]);
end

