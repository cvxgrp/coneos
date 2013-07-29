clear all

disp('------------------------------------------------------------')
disp('WARNING: this can take a very long time to run (due to CVX).')
disp('It may also crash/run out of memory.')
disp('Set run_cvx = false if you just want to use coneOS.')
disp('------------------------------------------------------------')

run ../coneOSdense/matlab/install_coneos_cvx.m

save_results = false;
run_cvx = false;
run_coneos = true;

randn('seed',0);rand('seed',0)

ns = [3000,10000,30000];
ms = ceil(ns/5);

time_pat_coneos = 'Time taken: (?<total>[\d\.]+)';
time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';
iter_pat_coneos = {'(?<iter>[\d]+)\|'};

for i = 1:length(ns)
    
    n=ns(i);
    m=ms(i);
    s=ceil(n/10);
    
    x_true=[randn(s,1);zeros(n-s,1)]; % true sparse signal
    x_true=x_true(randperm(n));
    A=randn(m,n);
    b = A*x_true + 0.1*randn(m,1); % measurements
    mu = 1;
    
    %%
    if run_coneos
        tic
        cvx_begin
        cvx_solver_settings('GEN_PLOTS',1) % only works for 'cvx_solver coneos_matlab'
        cvx_solver coneos
        variable x_c(n)
        minimize(0.5*sum_square(A*x_c - b) + mu*norm(x_c,1))
        output = evalc('cvx_end')
        toc
        
        coneos_direct.x{i} = x_c;
        coneos_direct.obj(i) = 0.5*sum_square(A*x_c - b) + mu*norm(x_c,1);
        timing = regexp(output, time_pat_coneos, 'names');
        coneos_direct.time{i} = str2num(timing.total);
        tmp = regexp(output, iter_pat_coneos, 'names');
        coneos_direct.iters{i} = str2num(tmp{1}(end).iter) + 1;
        coneos_direct.output{i} = output;


        if (save_results); save('data/lasso_coneos_direct', 'coneos_direct'); end

        %%
        
        tic
        cvx_begin
        cvx_solver_settings('USE_INDIRECT',1,'CG_MAX_ITS',2)
        cvx_solver_settings('GEN_PLOTS',1) % only works if 'cvx_solver coneos_matlab'
        cvx_solver coneos
        variable x_c(n)
        minimize(0.5*sum_square(A*x_c - b) + mu*norm(x_c,1))
        output = evalc('cvx_end')
        toc
        
        coneos_indirect.x{i} = x_c;
        coneos_indirect.obj(i) = 0.5*sum_square(A*x_c - b) + mu*norm(x_c,1);
        timing = regexp(output, time_pat_coneos, 'names');
        coneos_indirect.time{i} = str2num(timing.total);
        tmp = regexp(output, iter_pat_coneos, 'names');
        coneos_indirect.iters{i} = str2num(tmp{1}(end).iter) + 1;
        coneos_indirect.output{i} = output;
        
        if (save_results); save('data/lasso_coneos_indirect', 'coneos_indirect'); end

    end
    %%
    if run_cvx
        tic
        cvx_begin
        variable x_s(n)
        minimize(0.5*sum_square(A*x_s - b) + mu*norm(x_s,1))
        output = evalc('cvx_end')
        toc
        
        cvx.x{i} = x_s;
        cvx.obj(i) = 0.5*sum_square(A*x_s - b) + mu*norm(x_s,1);
        timing = regexp(output, time_pat_cvx, 'names');
        cvx.time{i} = str2num(timing.total);
        cvx.output{i} = output;
        
        if (save_results); save('data/lasso_cvx', 'cvx'); end;
    end
end
