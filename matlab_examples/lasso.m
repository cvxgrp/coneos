clear all
close all
randn('seed',0);rand('seed',0)
n=10000;
s=round(n/10);
m=1000;

x_true=[randn(s,1);zeros(n-s,1)]; % true sparse signal
x_true=x_true(randperm(n));
A=randn(m,n);
b = A*x_true + 0.1*randn(m,1); % measurements
%b = randn(m,1);
mu = 1;

cvx_begin
%cvx_solver_settings('NORMALIZE',0)
%cvx_solver_settings('USE_INDIRECT',1)
%cvx_solver_settings('ALPHA',1)
cvx_solver coneos
variable x_c(n)
minimize(0.5*sum_square(A*x_c - b) + mu*norm(x_c,1))
cvx_end

tic
x=zeros(n,1);y=x;theta=1;
t=1/max(eig(A*A')); % step-size
Ax=A*x;Ay=Ax;
fs=0.5*norm(Ax-b)^2 + mu*norm(x,1);
for k=1:3000
    if mod(k,500)==0
        k
    end
    x_old=x;
    Axold=Ax;
    
    temp=y-t*A'*(Ay-b);
    x = sign(temp).*max(abs(temp)-mu*t,0); % soft-thresholding
    Ax=A*x;
    
    theta_old=theta;
    theta = 0.5*(1+sqrt(1+4*theta^2));
    beta = (theta_old-1)/theta;
    
    y = x+beta*(x-x_old);
    Ay = (1+beta)*Ax - beta*Axold;
    fs=[fs;0.5*norm(Ax-b)^2 + mu*norm(x,1)];
    if (fs(end)>fs(end-1))
        y=x;
        theta=1;
    end
end
toc

%%
cvx_begin
variable x_s(n)
minimize(0.5*sum_square(A*x_s - b) + mu*norm(x_s,1))
cvx_end