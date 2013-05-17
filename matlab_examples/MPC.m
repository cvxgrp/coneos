clear all;
randn('state',0);
rand('state',0);

%% Defining the problem's parameters
% System's dimensions and control horizon
% n - states, m - controls
%n = 5; m = 2; T = 10; x_init = 5*randn(n,1);
n = 20; m = 5; T = 20; x_init = 5*randn(n,1);
%n = 100; m = 30; T = 50;  x_init = 5*randn(n,1);

A = randn(n);
A = A/max(abs(eig(A)));
B = randn(n,m);
B = 1.1*B / max(abs(svd(B)));
c = zeros(n,1);

mat = randn(n+m); mat = mat*mat';
mat = mat/norm(mat);
mat(1:n,n+1:end) = zeros(n,m);
mat(n+1:end,1:n) = zeros(m,n);
mat = sparse(mat);

q = zeros(n,1);
r = zeros(m,1);

umax = 1;
umin = -1;

cvx_begin
cvx_solver 'coneos'
%cvx_solver_settings('NORMALIZE',0)
cvx_solver_settings('ALPHA',1)
variables X_c(n,T+1) U_c(m,T+1)
obj=0;
for t=1:T+1
    obj = obj + 0.5 * [X_c(:,t);U_c(:,t)]'*mat*[X_c(:,t);U_c(:,t)]...
        + q'*X_c(:,t) + r'*U_c(:,t);
end
for t=1:T
    X_c(:,t+1) == A*X_c(:,t) + B*U_c(:,t) + c;
end
minimize (obj/T)
subject to
X_c(:,1) == x_init;
U_c <= umax;
U_c >= umin;
cvx_end


cvx_begin
variables X(n,T+1) U(m,T+1)
obj=0;
for t=1:T+1
    obj = obj + 0.5 * [X(:,t);U(:,t)]'*mat*[X(:,t);U(:,t)]...
        + q'*X(:,t) + r'*U(:,t);
end
for t=1:T
    X(:,t+1) == A*X(:,t) + B*U(:,t) + c;
end
minimize (obj/T)
subject to
X(:,1) == x_init;
U <= umax;
U >= umin;
cvx_end

