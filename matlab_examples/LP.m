clear all
close all
randn('seed',0);rand('seed',0)
n=1200;
m=3600;

A=randn(m,n);
b=rand(m,1);
c=rand(n,1);

tic
cvx_begin quiet
%cvx_solver_settings('USE_INDIRECT',1)
%cvx_solver_settings('NORMALIZE',0)
cvx_solver coneos
variable x_c(n)
dual variable z_c
minimize(c'*x_c)
z_c:A*x_c<=b
cvx_end
toc

tic
cvx_begin
%cvx_solver_settings('USE_INDIRECT',1)
%cvx_solver_settings('NORMALIZE',0)
cvx_solver pdos
variable x_p(n)
dual variable z_p
minimize(c'*x_p)
z_p:A*x_p<=b
cvx_end
toc

tic
cvx_begin
variable x(n)
minimize(c'*x)
A*x<=b
cvx_end
toc