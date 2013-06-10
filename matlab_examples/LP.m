clear all
close all
%%
randn('seed',0);rand('seed',0)
n=1800;
m=7200;

A=randn(m,n);
b=rand(m,1);
c=rand(n,1);

%%
tic
cvx_begin
%cvx_solver_settings('USE_INDIRECT',1)
cvx_solver_settings('NORMALIZE',1)
cvx_solver_settings('ALPHA',1.8)
cvx_solver_settings('RHO_X',1e-3)
cvx_solver coneos
variable x_c(n)
dual variable z_c
minimize(c'*x_c)
z_c:A*x_c<=b
cvx_end
toc

%%
cvx_begin
cvx_solver_settings('USE_INDIRECT',0)
cvx_solver_settings('CG_MAX_ITERS',1)
cvx_solver_settings('GEN_PLOTS',1)
cvx_solver_settings('RHOX',1e-3)
cvx_solver_settings('EPS',1e-6)

cvx_solver_settings('NORMALIZE',1)
cvx_solver_settings('ALPHA',1.0)
cvx_solver_settings('MAX_ITERS',1000)
cvx_solver_settings('PDOS_NORM',0)
%cvx_solver_settings('SIG',0.5*(1+sqrt(5)))
cvx_solver_settings('RELAX_X',0)
cvx_solver coneos_matlab
variable x_m(n)
dual variable z_m
minimize(c'*x_m)
z_m:A*x_m<=b
cvx_end

%%

tic
cvx_begin
cvx_solver pdos
variable x_p(n)
dual variable z_p
minimize(c'*x_p)
z_p:A*x_p<=b
cvx_end
toc
%%

tic
cvx_begin
variable x(n)
minimize(c'*x)
A*x<=b
cvx_end
toc