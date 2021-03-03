clear all;
close all;
syms f1 f2 f3
syms K11 K12 K13 K21 K22 K23 K33

K = [K11 K12 K13;
    K12 K22 K23;
    K13 K23 K33];

f=[f1 f2 f3].'

F=inv(eye(3)+hat(f))*(eye(3)-hat(f))

simplify(vee(F-F.') + 4*f/(1+f.'*f))


simplify( (trace(F)*eye(3)-F.')*(1+f.'*f) + 2*f*f.' + 2*hat(f) -2*eye(3) ) 


A = 1/2*(trace(F)*eye(3)-F.') ;

B = (eye(3) - hat(f) - f*f.')/(1+f.'*f); 

simplify(A-B)

%simplify(F*(1+f.'*f) - eye(3) +2*hat(f) -2*f*f.' + (f.'*f)*eye(3))

mu = 1/2*1/2*(trace(F)*eye(3)-F.')*K*vee(F-F.')

simplify(diff(mu,f1))