close all; clear all

run ../coneOSsparse/matlab/install_coneos_cvx.m
copyfile('../coneOSsparse/matlab/coneos_direct.m*','.');
copyfile('coneos_matlab/coneos_matlab.m','.');

cvx_on = true;
tests = dir('DIMACS/*.mat');
%params = struct('ALPHA',1.8, 'MAX_ITERS', 2000, 'VERBOSE', 1, 'NORMALIZE', 1, 'GEN_PLOTS',1);
params = struct('VERBOSE', 1);

time_pat_coneos = 'Time taken: (?<total>[\d\.]+)';
time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';
iter_pat_coneos = {'(?<iter>[\d]+)\|'};

%%

for i = 1:10%length(tests),
    
    %%%%%%%%%
    % NB: coneos solves the dual of these problems
    %%%%%%%%%
    
    clear A At b c K
    test_name = tests(i).name;
    f = ['DIMACS/' test_name];
    test_name = test_name(1:end-4); % strip .mat
    disp(['solving ' f]);
    load(f)
    
    [m1,n1] = size(b);
    [m2,n2] = size(c);
    if m1 == 1,
        b = b';
    end
    if m2 == 1,
        c = c';
    end
    if (exist('A','var'))
        data = struct('A', sparse(A'), 'b', full(c), 'c', -full(b));
    else
        data = struct('A', sparse(At), 'b', full(c), 'c', -full(b));
    end
    ccones = {'f','l','q','s'};
    for j = 1:length(ccones)
        if ~isfield(K,ccones{j})
            K.(ccones{j}) = [];
        end
    end
    
    cone = struct('f', 0, 'l', K.l, 'q', K.q', 's', K.s');
    
    [m1,n1] = size(cone.q);
    if m1 == 1,
        cone.q = cone.q';
    end
    [m1,n1] = size(cone.s);
    if m1 == 1,
        cone.s = cone.s';
    end
    
    if cvx_on
        % CVX:
        [m,n] = size(data.A);
        
        cvx_begin %quiet
        %cvx_solver coneos%_matlab
        variables xcvx(n) scvx(m)
        minimize(data.c'*xcvx)
        data.A*xcvx + scvx == data.b
        scvx(1:cone.f) == 0
        scvx(cone.f+1:cone.f+cone.l) >= 0
        idx=cone.f+cone.l;
        for kk =1:length(cone.q)
            norm(scvx(idx + 2: idx + cone.q(kk))) <= scvx(idx + 1)
            idx = idx + cone.q(kk);
        end
        for kk =1:length(cone.s)
            reshape(scvx(idx+1:idx + cone.s(kk)^2),cone.s(kk),cone.s(kk)) == semidefinite(cone.s(kk));
            idx = idx + cone.s(kk)^2;
        end
        cvx_end
        
        cvx_objval.(test_name) =  cvx_optval;
    end
    
    tic
    [output, x_m, z_m] = evalc('coneos_matlab(data,cone,params);');
    output
    toc
    
    coneos_NNZA.(test_name) = nnz(data.A);
    coneos_x.(test_name) = x_m;
    timing = regexp(output, time_pat_coneos, 'names');
    coneos_times.(test_name) = str2num(timing.total);
    tmp = regexp(output,iter_pat_coneos, 'names');
    coneos_iters.(test_name) = str2num(tmp{1}(end).iter) + 1;
    coneos_objval.(test_name) = data.c'*x_m;
    coneos_output.(test_name) = output;
    save data/dimacs_run_data
    
end