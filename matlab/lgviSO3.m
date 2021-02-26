%Lie group variational integrator
%for the attitude dynamics of a 3D Pendulum on SO(3)
clear all;
close all;

%Pendulum Properties 
m=1;
g=0;
e3=[0 0 1].';

J=diag([1 2.8 2]);
Jd=trace(J)/2*eye(3)-J;%nonstandard inertia matrix
rho=e3;

%Simulation time
T=20;
h=0.01;
N=T/h+1;
t=linspace(0,T,N);

%Initial conditions
R(:,:,1)=diag([-1,1,-1]);%rotation matrix
%R(:,:,1)=[0 0 1;
    %0 1 0;
    %-1 0 0];
W(:,1)=[0.5 -0.5 0.4];%angular velocity


%Lie group variaitonal integrator
M(:,1)=m*g*cross(rho,R(:,:,1)'*e3);%moment due to gravity

for k=1:N-1
 
    %Newton iteration 
    gk=h*J*W(:,k)+h^2/2*M(:,k);
    Fk = findF(gk, J, h*W(:,k));

    %update R,M,w
    R(:,:,k+1)=R(:,:,k)*Fk;
    M(:,k+1)=m*g*cross(rho,R(:,:,k+1)'*e3);
    W(:,k+1)=J\(Fk'*(J*W(:,k)+h/2*M(:,k))+h/2*M(:,k+1));
end

%Total energy 
for k=1:N
    E(k)=1/2*W(:,k)'*J*W(:,k)-m*g*rho'*R(:,:,k)'*e3;
end


plot(t,E);
ylabel('$$E$$','interpreter','latex');
xlabel('$$t$$','interpreter','latex');

save lgviSO3

