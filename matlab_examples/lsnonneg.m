clear all
close all
%%
randn('seed',124);rand('seed',123);
m = 1000;
n = 1000;
A = randn(m,n);
x_true = pos(randn(n,1));
%b = A*x_true + 0.1*randn(m,1);
b=randn(m,1);
mu = 1;

%%
tic
cvx_begin
cvx_solver coneos_matlab
cvx_solver_settings('USE_INDIRECT',1)
cvx_solver_settings('CG_MAX_ITERS',1)
cvx_solver_settings('GEN_PLOTS',1)
cvx_solver_settings('RHOX',1e-3)
%cvx_solver_settings('CG_TOL',1e-4)
cvx_solver_settings('NORMALIZE',1)
cvx_solver_settings('ALPHA',1.8)
cvx_solver_settings('MAX_ITERS',2000)
%cvx_solver_settings('SIG',0.5*(1+sqrt(5)))
cvx_solver_settings('RELAX_X',0)
cvx_solver_settings('PDOS_NORM',0)
variable x_m(n)
dual variable y_m
minimize(sum_square(A*x_m - b) + mu*sum_square(x_m))
y_m:x_m>=0
cvx_end
toc
%%

tic
cvx_begin
cvx_solver coneos
%cvx_solver_settings('NORMALIZE',0)
cvx_solver_settings('USE_INDIRECT',0)
cvx_solver_settings('CG_MAX_ITS',2)
cvx_solver_settings('CG_TOL',1e-12)
%cvx_solver_settings('NORMALIZE',1)
%cvx_solver_settings('ALPHA',1.8)
%cvx_solver_settings('RHO_X',1e-3)
%cvx_solver_settings('EPS',1e-8)
cvx_solver_settings('MAX_ITERS',5000)
variable x_c(n)
dual variable y_c
minimize(sum_square(A*x_c - b) + mu*sum_square(x_c))
y_c:x_c>=0
cvx_end
toc

%%
tic
cvx_begin
cvx_solver pdos
variable x_p(n)
dual variable y_p
minimize(sum_square(A*x_p - b) + mu*sum_square(x_p))
y_p:x_p>=0
cvx_end
toc

%%
tic
cvx_begin
variable x_cvx(n)
dual variable y_cvx
minimize(sum_square(A*x_cvx - b) + mu*sum_square(x_cvx))
y_cvx:x_cvx>=0
cvx_end
toc
%%

% solves: min. norm(Ax - b)^2 + mu*norm(x)^2, s.t. x>=0
tic
rho=1;
[L,D,P] = ldl(A'*A+(rho+mu)*eye(n),'vector');
z = zeros(n,1);lam = z;
Ab = A'*b;
iter = 1;
while (iter < 10000)
    if (mod(iter,100)==0)
        iter
    end
    tmp = Ab+rho*(z-lam);
    x = (L'\(D\(L\tmp(P))));x(P)=x;
    z = pos(x+lam);
    lam = lam + x - z;
    if norm(x-z) < 1e-4; break; end
    iter = iter+1;
end
toc