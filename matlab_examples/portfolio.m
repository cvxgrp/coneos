clear all

disp('------------------------------------------------------------')
disp('WARNING: this can take a very long time to run.')
disp('It may also crash/run out of memory.')
disp('Set run_sdpt3 = false if you just want to run coneOS.')
disp('------------------------------------------------------------')

run ../coneOSsparse/matlab/install_coneos_cvx.m

save_results = false;
run_sdpt3 = false;
run_coneos = true;

ns = [5000, 50000, 100000];
ms = [50, 500, 1000];

time_pat_coneos = 'Time taken: (?<total>[\d\.]+)';
time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';
iter_pat_coneos = {'(?<iter>[\d]+)\|'};

for i = 1:length(ns)
    seedstr = sprintf('coneos_portfolio_ex_%i',i);
    randn('seed',sum(seedstr));rand('seed',sum(seedstr))
    
    n = ns(i);
    m = ms(i);
    
    mu = exp(randn(n,1));
    D = sqrt(2*rand(n,1));
    F = randn(n,m);
    gamma = 10;
    
    portfolio_prob.F{i} = F;
    portfolio_prob.D{i} = D;
    portfolio_prob.mu{i} = mu;
    portfolio_prob.n{i} = n;
    portfolio_prob.m{i} = m;
    portfolio_prob.gamma{i} = gamma;
    
    if (save_results); save('data/portfolio_prob', 'portfolio_prob','-v7.3'); end
    
    %%
    if run_coneos
        
        tic
        cvx_begin
        cvx_solver coneos
        cvx_solver_settings('GEN_PLOTS',1) % only works if 'cvx_solver coneos_matlab'
        variable x(n)
        maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
        sum(x) == 1
        x >= 0
        output = evalc('cvx_end')
        toc
        
        coneos_direct.x{i} = x;
        coneos_direct.x_viol{i} = min(x);
        coneos_direct.budget_viol{i} = abs(1-sum(x));
        coneos_direct.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
        timing = regexp(output, time_pat_coneos, 'names');
        coneos_direct.time{i} = str2num(timing.total);
        tmp = regexp(output, iter_pat_coneos, 'names');
        coneos_direct.iters{i} = str2num(tmp{1}(end).iter) + 1;
        coneos_direct.output{i} = output;
        
        if (save_results); save('data/portfolio_coneos_direct', 'coneos_direct'); end
        
        %%
        
        tic
        cvx_begin
        cvx_solver coneos
        cvx_solver_settings('USE_INDIRECT',1,'CG_MAX_ITS',2)
        cvx_solver_settings('GEN_PLOTS',1) % only works if 'cvx_solver coneos_matlab'
        variable x(n)
        maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
        sum(x) == 1
        x >= 0
        output = evalc('cvx_end')
        toc
        
        coneos_indirect.x{i} = x;
        coneos_indirect.x_viol{i} = min(x);
        coneos_indirect.budget_viol{i} = abs(1-sum(x));
        coneos_indirect.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
        timing = regexp(output, time_pat_coneos, 'names');
        coneos_indirect.time{i} = str2num(timing.total);
        tmp = regexp(output, iter_pat_coneos, 'names');
        coneos_indirect.iters{i} = str2num(tmp{1}(end).iter) + 1;
        coneos_indirect.output{i} = output;
        
        
        if (save_results); save('data/portfolio_coneos_indirect', 'coneos_indirect'); end
        
    end
    %%
    if run_sdpt3
        try
            tic
            cvx_begin
            cvx_solver sdpt3
            variable x(n)
            maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
            sum(x) == 1
            x >= 0
            output = evalc('cvx_end')
            toc
            
            cvx.x{i} = x;
            cvx.x_viol{i} = min(x);
            cvx.budget_viol{i} = abs(1-sum(x));
            cvx.obj(i) = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
            timing = regexp(output, time_pat_cvx, 'names');
            cvx.time{i} = str2num(timing.total);
            cvx.output{i} = output;
            cvx.err{i} = 0;
            
        catch err
            cvx.time{i} = toc;
            cvx.err{i} = err;
        end
        
        if (save_results); save('data/portfolio_cvx', 'cvx'); end
        
    end
end