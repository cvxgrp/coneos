clear all
close all
%%
randn('seed',0);rand('seed',0);
n=8000;
m=10000;

A=randn(m,n);
b=rand(m,1);
p=n*2;
R=randn(p,n);
c=rand(p,1);

%%
tic
cvx_begin
cvx_solver 'coneos'
variable x_c(n)
cvx_solver_settings('MAX_ITERS',2000)
%cvx_solver_settings('USE_INDIRECT',1)
cvx_solver_settings('NORMALIZE',1)
cvx_solver_settings('CG_MAX_ITS',1)
cvx_solver_settings('CG_TOL',1e-12)
dual variable z_c
minimize(0.5*sum_square(R*x_c - c))
z_c:A*x_c<=b
cvx_end
toc

%%

tic
cvx_begin
%cvx_solver_settings('USE_INDIRECT',1)
cvx_solver_settings('USE_INDIRECT',1)
cvx_solver_settings('CG_MAX_ITERS',2)
cvx_solver_settings('CG_TOL',1e-12)
cvx_solver_settings('GEN_PLOTS',1)
cvx_solver_settings('RHOX',1e-3)
cvx_solver_settings('NORMALIZE',1)
cvx_solver_settings('ALPHA',1.8)
cvx_solver_settings('MAX_ITERS',2000)
cvx_solver_settings('EPS',1e-5)
%cvx_solver_settings('SIG',0.5*(1+sqrt(5)))
cvx_solver_settings('RELAX_X',0)
cvx_solver 'coneos_matlab'
variable x_c(n)
dual variable z_c
minimize(0.5*sum_square(R*x_c - c))
z_c:A*x_c<=b
cvx_end
toc
%%
tic
cvx_begin
cvx_solver 'pdos'
variable x_p(n)
dual variable z_p
minimize(0.5*sum_square(R*x_p - c))
z_p:A*x_p<=b
cvx_end
toc

%%
tic
cvx_begin
cvx_solver 'sdpt3'
variable x_cvx(n)
dual variable z_cvx
minimize(0.5*sum_square(R*x_cvx - c))
z_cvx:A*x_cvx<=b
cvx_end
toc