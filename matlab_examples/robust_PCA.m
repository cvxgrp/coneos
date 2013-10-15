clear all

disp('------------------------------------------------------------')
disp('WARNING: this can take a very long time to run.')
disp('It may also crash/run out of memory.')
disp('Set run_sdpt3 = false if you just want to run coneOS.')
disp('------------------------------------------------------------')

run ../coneOSsparse/matlab/install_coneos_cvx.m

save_results = true;
run_sdpt3 = false;
run_coneos = true;

ns = [100,500,1000];
ms = ns; % square matrices, but doesn't have to be

time_pat_coneos = 'Time taken: (?<total>[\d\.]+)';
time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';
iter_pat_coneos = {'(?<iter>[\d]+)\|'};

for i = 1:length(ns)
    seedstr = sprintf('coneos_rpca_ex_%i',i);
    randn('seed',sum(seedstr));rand('seed',sum(seedstr))

    n = ns(i);
    m = ms(i);
    r = 10; % rank
    
    L1 = randn(m,r);
    L2 = randn(r,n);
    L = L1*L2;
    S = sprandn(m,n,0.05);
    M = L + S;
    kap = sum(norms(S,1));
    
    rpca_prob.L{i} = L;
    rpca_prob.S{i} = S;
    rpca_prob.M{i} = M;
    rpca_prob.kap{i} = kap;
    rpca_prob.r{i} = r;
    rpca_prob.n{i} = n;
    rpca_prob.m{i} = m;
    if (save_results); save('data/rpca_prob', 'rpca_prob','-v7.3'); end
    
    %%
    if run_coneos
        
        tic
        cvx_begin
        cvx_solver coneos
        cvx_solver_settings('GEN_PLOTS',1) % only works if 'cvx_solver coneos_matlab'
        variables Lc(m,n) Sc(m,n)
        dual variable Yc
        minimize(norm_nuc(Lc))
        sum(norms(Sc,1)) <= kap
        Yc:Lc + Sc == M;
        output = evalc('cvx_end')
        toc
        
        coneos_direct.L{i} = Lc;
        coneos_direct.obj(i) = norm_nuc(Lc);
        timing = regexp(output, time_pat_coneos, 'names');
        coneos_direct.time{i} = str2num(timing.total);
        tmp = regexp(output, iter_pat_coneos, 'names');
        coneos_direct.iters{i} = str2num(tmp{1}(end).iter) + 1;
        coneos_direct.output{i} = output;
        
        
        if (save_results); save('data/rpca_coneos_direct', 'coneos_direct'); end
        
        %%
        
        tic
        cvx_begin
        cvx_solver coneos
        cvx_solver_settings('USE_INDIRECT',1,'CG_MAX_ITS',5)
        cvx_solver_settings('GEN_PLOTS',1) % only works if 'cvx_solver coneos_matlab'
        variables Lc(m,n) Sc(m,n)
        dual variable Yc
        minimize(norm_nuc(Lc))
        sum(norms(Sc,1)) <= kap
        Yc:Lc + Sc == M;
        output = evalc('cvx_end')
        toc
        
        coneos_indirect.L{i} = Lc;
        coneos_indirect.obj(i) = norm_nuc(Lc);
        timing = regexp(output, time_pat_coneos, 'names');
        coneos_indirect.time{i} = str2num(timing.total);
        tmp = regexp(output, iter_pat_coneos, 'names');
        coneos_indirect.iters{i} = str2num(tmp{1}(end).iter) + 1;
        coneos_indirect.output{i} = output;
        
        if (save_results); save('data/rpca_coneos_indirect', 'coneos_indirect'); end
        
    end
    %%
    if run_sdpt3
        try
            tic
            cvx_begin
            cvx_solver sdpt3
            variables Lt(m,n) St(m,n)
            dual variable Yt
            minimize(norm_nuc(Lt))
            sum(norms(St,1)) <= kap
            Yt:Lt + St == M;
            output = evalc('cvx_end')
            toc
            
            cvx.L{i} = Lt;
            cvx.obj(i) = norm_nuc(Lt);
            timing = regexp(output, time_pat_cvx, 'names');
            cvx.time{i} = str2num(timing.total);
            cvx.output{i} = output;
            cvx.err{i} = 0;
            
        catch err
            cvx.time{i} = toc;
            cvx.err{i} = err;
        end
        
        if (save_results); save('data/rpca_cvx', 'cvx'); end
        
    end
end
