clear all
close all
%randn('seed',0);rand('seed',0)
n=4000;
s=round(n/10);
m=300;

x_true=[randn(s,1);zeros(n-s,1)]; % true sparse signal
x_true=x_true(randperm(n));
A=randn(m,n);
b = A*x_true + 0.1*randn(m,1); % measurements

mu = 1;

cvx_begin
cvx_solver coneos
variable x_c(n)
minimize(sum_square(A*x_c - b) + mu*norm(x_c,1))
cvx_end

cvx_begin
variable x(n)
minimize(sum_square(A*x - b) + mu*norm(x,1))
cvx_end