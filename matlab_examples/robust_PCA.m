n = 500;
r = 10;
L = randn(r,n);
L = L'*L;
S = sprandn(n,n,0.05);
M = L + S;
mu = 0.1;

cvx_begin
cvx_solver coneos
%cvx_solver_settings('NORMALIZE',0)
%cvx_precision high
variables Lc(n,n) Sc(n,n)
minimize(norm_nuc(Lc) + mu*sum(norms(Sc,1)))
Lc + Sc == M;
cvx_end


cvx_begin
variables Lt(n,n) St(n,n)
minimize(norm_nuc(Lt) + mu*sum(norms(St,1)))
Lt + St == M;
cvx_end