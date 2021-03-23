clear all;
close all;

syms R11 R12 R13 R21 R22 R23 R31 R32 R33 x1 x2 x3 K11 K22 K13 K23 P1 P2 P3 e1 e2 e3
K = [K11 0 K13;
    0 K11 K23;
    0 0 1];
R=[R11 R12 R13;
    R21 R22 R23;
    R31 R32 R33];
P = [P1 P2 P3 1].';
x= [x1 x2 x3].';
eta =[e1 e2 e3].';

p = K*[R*hat(eta), x]*P;
uv = p(1:2)/p(3);

gradx1=simplify(diff(uv,x1));
gradx2=simplify(diff(uv,x2));
gradx3=simplify(diff(uv,x3));

gradR1 = simplify(subs(diff(uv,e1),[e1,e2,e3],[0,0,0]));
gradR2 = simplify(subs(diff(uv,e2),[e1,e2,e3],[0,0,0]));
gradR3 = simplify(subs(diff(uv,e3),[e1,e2,e3],[0,0,0]));

