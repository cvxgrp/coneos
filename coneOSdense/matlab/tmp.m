clear all
close all
%randn('seed',0);rand('seed',0
n=40;
m=120;

A=randn(m,n);
b=rand(m,1);
c=rand(n,1);


cone = struct('f',0,'l',m,'k_soc',0,'q',[],'ssize',0,'s',[]); 
data = struct('A',full(A),'b',b,'c',c);
params = struct('ALPHA',1.8,'UNDET_TOL',1e-4,'EPS_ABS',1e-4,'EPS_REL',1e-4,'MAX_ITERS',2000,...
            'CG_MAX_ITS',20,'CG_TOL',1e-3,'GEN_PLOTS',1,'VERBOSE',1,'NORMALIZE',0);
        
[x_d,z_d,status_d] = coneOS_direct(data,cone,params);
%[x_d,z_d,status_d] = coneOS_indirect(data,cone,params);


cvx_begin
cvx_solver coneos
variable x(n)
dual variable z
minimize(c'*x)
z:A*x<=b
cvx_end

cvx_begin
variable x(n)
minimize(c'*x)
A*x<=b
cvx_end