clear all;close all
randn('seed',0);rand('seed',0);
m = 80;
n = 50;
A = randn(m,n);
b = randn(m,1);
mu = 1;

 
%% CONEOS
tic
cvx_begin
cvx_solver coneos
variable x_c(n)
minimize(sum_square(A*x_c - b) + mu*sum_square(x_c))
cvx_end
toc

%% PDOS
tic
cvx_begin
cvx_solver pdos
variable x_p(n)
minimize(sum_square(A*x_p - b) + mu*sum_square(x_p))
cvx_end
toc

%% SDPT3
tic
cvx_begin
cvx_solver sdpt3
variable x_s(n)
minimize(sum_square(A*x_s - b) + mu*sum_square(x_s))
cvx_end
toc 

