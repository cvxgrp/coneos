clear all

disp('------------------------------------------------------------')
disp('WARNING: this can take a very long time to run (due to CVX).')
disp('It may also crash/run out of memory.')
disp('Set run_cvx = false if you just want to use coneOS.')
disp('------------------------------------------------------------')

run ../coneOSsparse/matlab/install_coneos_cvx.m

save_results = true;
run_cvx = false;
run_coneos = true;

randn('seed',0);rand('seed',0)

ns = [50,150,300];
ms = ns; % square matrices, but doesn't have to be

time_pat_coneos = 'Time taken: (?<total>[\d\.]+)';
time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';
iter_pat_coneos = {'(?<iter>[\d]+)\|'};

for i = 1:length(ns)
    
    n = ns(i);
    m = ms(i);
    r = 10; % rank
    
    L1 = randn(m,r);
    L2 = randn(r,n);
    L = L1*L2;
    S = sprandn(m,n,0.05);
    M = L + S;
    kap = sum(norms(S,1));

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
        cvx_solver_settings('USE_INDIRECT',1,'CG_MAX_ITS',2)
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
    if run_cvx
        tic
        cvx_begin
        variables Lt(m,n) St(m,n)
        dual variable Yt
        minimize(norm_nuc(Lt))
        sum(norms(St,1)) <= kap
        Yt:Lt + St == M;
        output = evalc('cvx_end')
        toc
        
        cvx.L{i} = Lt;
        cvx.obj(i) = norm_nuc(Lc);
        timing = regexp(output, time_pat_cvx, 'names');
        cvx.time{i} = str2num(timing.total);
        cvx.output{i} = output;
        
        if (save_results); save('data/rpca_cvx', 'cvx'); end
    end
end
