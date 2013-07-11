clear all; close all
n = 50; m=n;
randn('seed',0);rand('seed',0)
r = 10;
L1 = randn(m,r);
L2 = randn(r,n);
L = L1*L2;
S = sprandn(m,n,0.05);
M = L + S;
mu = 0.1;
kap = sum(norms(S,1));

%%
for i=1:50
cvx_begin
cvx_solver coneos
cvx_solver_settings('USE_INDIRECT',1)
cvx_solver_settings('CG_MAX_ITS',1)
cvx_solver_settings('CG_TOL',1e-12)
cvx_solver_settings('GEN_PLOTS',1)
cvx_solver_settings('EPS',1e-6)
cvx_solver_settings('MAX_ITERS',5000)
variables Lc(m,n) Sc(m,n)
dual variable Yc
minimize(norm_nuc(Lc))
sum(norms(Sc,1)) <= kap
Yc:Lc + Sc == M;
cvx_end
end

%%
tic
cvx_begin
variables Lt(m,n) St(m,n)
dual variable Yt
minimize(norm_nuc(Lt))
sum(norms(St,1)) <= kap
Yt:Lt + St == M;
cvx_end
toc