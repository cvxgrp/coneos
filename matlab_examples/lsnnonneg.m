clear all
close all

randn('seed',0);rand('seed',0);
m = 8000;
n = 5000;
A = randn(m,n);
x_true = pos(randn(n,1));
b = A*x_true + 0.1*randn(m,1);
%b=randn(m,1);
mu = 1;

tic
cvx_begin
cvx_solver coneos
variable x_c(n)
dual variable y_c
minimize(sum_square(A*x_c - b) + mu*sum_square(x_c))
y_c:x_c>=0
cvx_end
toc

tic
cvx_begin
cvx_solver pdos
variable x_p(n)
dual variable y_p
minimize(sum_square(A*x_p - b) + mu*sum_square(x_p))
y_p:x_p>=0
cvx_end
toc

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
while 1
    tmp = Ab+rho*(z-lam);
    x = (L'\(D\(L\tmp(P))));x(P)=x;
    z = pos(x+lam);
    lam = lam + x - z;
    if norm(x-z) < 1e-4; break; end
end
toc