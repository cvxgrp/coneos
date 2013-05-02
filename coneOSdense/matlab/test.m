clear all
close all
%randn('seed',0);rand('seed',0
n=30;
m=5;
mt=m^2;

m2 = 4;
mt2= m2^2;

A=randn(mt,n);
b=rand(m,m);
b=(b+b')/2;
b=reshape(b,mt,1);
c=rand(n,1);

A2=randn(mt2,n);
b2=rand(m2,m2);
b2=(b2+b2')/2;
b2=reshape(b2,mt2,1);


tic
cvx_begin
cvx_solver coneos
variable x(n) 
variable S(m,m) symmetric
variable S2(m2,m2) symmetric
dual variables yy yy2
minimize(c'*x)
yy:S==semidefinite(m);
yy2:S2==semidefinite(m2);
x(1:3)>=0
A*x + S(:) == b
A2*x + S2(:) == b2

cvx_end
toc

%%
n = 200;
m = 20;
A = randn(m,n);
b = randn(m,1);
mu = 1;

tic
cvx_begin
cvx_solver coneos
variable x(n)
minimize(sum_square(A*x-b)+mu*norm(x,1))
cvx_end
toc

tic
cvx_begin
cvx_solver
variable x_sed(n)
minimize(sum_square(A*x_sed-b)+mu*norm(x_sed,1))
cvx_end
toc


%%

clear all
close all
randn('seed',1); rand('seed',1);

n = 12; m = 4;l=n+m+1;
A = randn(n,n); A=A/max(abs(eig(A)));
B = randn(n,m); B=B/max(svd(B));

Q = randn(n,n);Q=10*(Q*Q')/n;
R = randn(m,m);R=(R*R')/n;

u_max = 1;

W = 0.1*randn(n,n);W=(W*W');
Wh = zeros(n+1);Wh(1:n,1:n) = W;
X_sig=0.2*eye(n);
x_mean = 10*randn(n,1);
Xh = [X_sig+x_mean*x_mean' x_mean;x_mean' 1];
gamma = 0.99;

F = sparse([A B zeros(n,1);zeros(1,n) zeros(1,m) 1]);
G = sparse([eye(n) zeros(n,m) zeros(n,1);zeros(1,n) zeros(1,m) 1]);
L = sparse([Q zeros(n,m) zeros(n,1);zeros(m,n) R zeros(m,1);zeros(1,n) zeros(1,m) 0]);
el = sparse(zeros(l,1));el(n+m+1)=1;

%%

T=1;
cvx_begin
cvx_solver coneos
variable P(n+1,n+1,T+1) symmetric
variables lambda(m,T)
dual variable X
for t=1:T
    %P(1:n,1:n,t) == semidefinite(n)   
    M=[zeros(n) zeros(n,m) zeros(n,1);zeros(m,n) diag(lambda(:,t)) zeros(m,1);zeros(1,n) zeros(1,m) -u_max^2*sum(lambda(:,t))];
    (L+gamma*F'*P(:,:,t+1)*F - G'*P(:,:,t)*G + M + gamma*trace(P(:,:,t+1)*Wh)*el*el')==semidefinite(l);
end
X:P(1:n,1:n,1) == semidefinite(n) 
lambda>=0
P(:,:,1)==P(:,:,T+1)
maximize(trace(P(:,:,1)*Xh))
cvx_end
