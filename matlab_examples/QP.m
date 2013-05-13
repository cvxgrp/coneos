clear all
close all
n=500;
m=1200;

A=randn(m,n);
b=rand(m,1);
p=n*2;
R=randn(p,n);
c=rand(p,1);

%%

cvx_begin
cvx_solver 'coneos'
variable x_cvx(n)
dual variable z_cvx
minimize(0.5*sum_square(R*x_cvx - c))
z_cvx:A*x_cvx<=b
cvx_end


cvx_begin
cvx_solver 'sdpt3'
variable x_cvx(n)
dual variable z_cvx
minimize(0.5*sum_square(R*x_cvx - c))
z_cvx:A*x_cvx<=b
cvx_end
