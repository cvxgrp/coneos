clear all

run ../coneOSdense/matlab/install_coneos_cvx.m

randn('seed',0);rand('seed',0)

ns = [3000,10000,30000];
ms = ceil(ns/5);

for i = 1:1%length(ns)
    
    n=ns(i);
    m=ms(i);
    s=ceil(n/10);
    
    x_true=[randn(s,1);zeros(n-s,1)]; % true sparse signal
    x_true=x_true(randperm(n));
    A=randn(m,n);
    b = A*x_true + 0.1*randn(m,1); % measurements
    mu = 1;
    
    %%
    
    tic
    cvx_begin
    cvx_solver_settings('GEN_PLOTS',1) % only works for 'cvx_solver coneos_matlab'
    cvx_solver coneos
    variable x_c(n)
    minimize(0.5*sum_square(A*x_c - b) + mu*norm(x_c,1))
    cvx_end
    toc
    
    xcd{i} = x_c;
    objcd(i) = 0.5*sum_square(A*x_c - b) + mu*norm(x_c,1);
    
    %%
    
    tic
    cvx_begin
    cvx_solver_settings('USE_INDIRECT',1,'CG_MAX_ITS',1)
    cvx_solver_settings('GEN_PLOTS',1) % only works if 'cvx_solver coneos_matlab'
    cvx_solver coneos
    variable x_c(n)
    minimize(0.5*sum_square(A*x_c - b) + mu*norm(x_c,1))
    cvx_end
    toc
    
    xci{i} = x_c;
    objci(i) = 0.5*sum_square(A*x_c - b) + mu*norm(x_c,1);

    %%
    
    tic
    cvx_begin
    variable x_s(n)
    minimize(0.5*sum_square(A*x_s - b) + mu*norm(x_s,1))
    cvx_end
    toc
    
    xs{i} = x_s;
    objs(i) = 0.5*sum_square(A*x_s - b) + mu*norm(x_s,1);
end