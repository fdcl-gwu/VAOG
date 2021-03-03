function FCay=findF(g,J,varargin)

if isempty(varargin) 
    f = zeros(3,1);
else
    f = varargin{1};
end

delf=1;
while delf > 1e-15
    G=g+(hat(g)-2*J)*f+(g'*f)*f;
    nabG=(hat(g)-2*J)+g'*f*eye(3)+f*g';
    f_new=f-inv(nabG)*G;
    delf=norm(f_new-f);
    f=f_new;
end


%FCay=(eye(3)+skew(f))*inv(eye(3)-skew(f));
FCay=[ 1+f(1)^2-f(3)^2-f(2)^2,     -2*f(3)+2*f(2)*f(1),      2*f(3)*f(1)+2*f(2);
    2*f(3)+2*f(2)*f(1), -f(3)^2+1+f(2)^2-f(1)^2,      2*f(3)*f(2)-2*f(1);
    -2*f(2)+2*f(3)*f(1),      2*f(3)*f(2)+2*f(1), -f(2)^2-f(1)^2+1+f(3)^2]/(1+f'*f);
