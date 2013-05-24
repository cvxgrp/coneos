n = 50; m=30;
randn('seed',0);rand('seed',0)
r = 10;
L1 = randn(m,r);
L2 = randn(r,n);
L = L1*L2;
S = sprandn(m,n,0.05);
M = L + S;
mu = 0.1;

cvx_begin
cvx_solver coneos
cvx_solver_settings('USE_INDIRECT',1)
%cvx_solver_settings('CG_TOL',1e-8)
%cvx_solver_settings('CG_MAX_ITS',100)
variables Lc(m,n) Sc(m,n)
dual variable Yc
minimize(norm_nuc(Lc) + mu*sum(norms(Sc,1)))
Yc:Lc + Sc == M;
cvx_end


cvx_begin
variables Lt(m,n) St(m,n)
minimize(norm_nuc(Lt) + mu*sum(norms(St,1)))
Lt + St == M;
cvx_end