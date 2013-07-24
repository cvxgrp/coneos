clear all

run ../coneOSsparse/matlab/install_coneos_cvx.m

randn('seed',0);rand('seed',0)

ns = [50,150,300];
ms = ns; % square matrices, but doesn't have to be

for i=1:1
    
    n = ns(i);
    m = ms(i);
    r = 10;
    
    L1 = randn(m,r);
    L2 = randn(r,n);
    L = L1*L2;
    S = sprandn(m,n,0.05);
    M = L + S;
    kap = sum(norms(S,1));
    
    %%
    tic
    cvx_begin
    cvx_solver coneos
    variables Lc(m,n) Sc(m,n)
    dual variable Yc
    minimize(norm_nuc(Lc))
    sum(norms(Sc,1)) <= kap
    Yc:Lc + Sc == M;
    cvx_end
    toc
    
    Lcd{i} = Lc;
    Scd{i} = Sc;
    objcd(i) = norm_nuc(Lc);

    
    %%
    tic
    cvx_begin
    cvx_solver coneos
    cvx_solver_settings('USE_INDIRECT',1,'CG_MAX_ITS',1)
    variables Lc(m,n) Sc(m,n)
    dual variable Yc
    minimize(norm_nuc(Lc))
    sum(norms(Sc,1)) <= kap
    Yc:Lc + Sc == M;
    cvx_end
    toc
    
    Lci{i} = Lc;
    Sci{i} = Sc;
    objci(i) = norm_nuc(Lc);

    
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
    
    Ls{i} = Lt;
    Ss{i} = St;
    objs(i) = norm_nuc(Lt);

end