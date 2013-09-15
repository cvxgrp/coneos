close all; clear all
cd 'DIMACS'

run ../../coneOSsparse/matlab/install_coneos_cvx.m
copyfile('../../coneOSsparse/matlab/coneos_direct.m*','.');

cvx_on = false;
coneos_on = true;
save_data = true;
tests = dir('*.mat');
params = struct('VERBOSE', 1, 'EPS_ABS', 1e-5, 'MAX_ITERS', 5000);
for i=1:length(tests)
    szes(i) = tests(i).bytes;
end
[szes, solve_order] = sort(szes);

time_pat_coneos = 'Time taken: (?<total>[\d\.]+)';
time_pat_cvx = 'Total CPU time \(secs\)\s*=\s*(?<total>[\d\.]+)';
iter_pat_coneos = {'(?<iter>[\d]+)\|'};

%%
N = length(tests);
for ii = 1:N
    i = solve_order(ii); %% solve in increasing order of size
    
    clear A At b c K
    test_name = tests(i).name;
    f = ['DIMACS/' test_name];
    test_name = test_name(1:end-4); % strip .mat
    fprintf('running test %i out of %i : %s\n', ii, N, test_name);

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
    if isempty(K.f)
        K.f = 0;
    end
    if isempty(K.l)
        K.l = 0;
    end
    if (K.q == 0)
        K.q = [];
    end
    if (K.s == 0)
        K.s = [];
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
    if (data.b == 0); data.b = zeros(size(data.A,1),1); end
    if (data.c == 0); data.c = zeros(size(data.A,2),1); end
    
    
    if cvx_on
        % CVX:
        [m,n] = size(data.A);
        
        cvx_begin %quiet
        cvx_solver coneos%_matlab
        variables xcvx(n) scvx(m)
        dual variable zcvx
        minimize(data.c'*xcvx)
        zcvx: data.A*xcvx + scvx == data.b
        scvx(1:cone.f) == 0
        scvx(cone.f+1:cone.f+cone.l) >= 0
        idx=cone.f+cone.l;
        if isempty(idx) idx = 0; end
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
    
    if coneos_on
        if (isfield(K,'r') && K.r ~= 0)
            coneos_error(test_name) = 'rotated lorentz cones not currently supported';
            coneos_x = nan;
            coneos_objval = nan;
        else
            
            tic
            [output, x_m, z_m] = evalc('coneos_direct(data,cone,params);');
            output
            data.b'*z_m %% because coneos targets dual of DIMACs formulations
            toc
            
            coneos_NNZA.(test_name) = nnz(data.A);
            coneos_x.(test_name) = x_m;
            coneos_z.(test_name) = z_m;
            timing = regexp(output, time_pat_coneos, 'names');
            coneos_times.(test_name) = str2num(timing.total);
            tmp = regexp(output,iter_pat_coneos, 'names');
            coneos_iters.(test_name) = str2num(tmp{1}(end).iter) + 1;
            coneos_dobjval.(test_name) = data.b'*z_m;
            coneos_pobjval.(test_name) = data.c'*x_m;
            coneos_output.(test_name) = output;
        end
    end
    
    if save_data; save ../data/dimacs_run_data; end
    
end
delete 'coneos_direct.m*'
cd ..