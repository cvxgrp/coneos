close all; clear all

cvx_on = true;

tests = dir('DIMACS/*.mat');

params = struct('ALPHA',1.8,'EPS_ABS', 1e-32, 'EPS_REL', 1e-32, 'CG_MAX_ITS', 1, 'CG_TOL', 1e-8, 'MAX_ITERS',5000, 'VERBOSE', 1, 'NORMALIZE', 1, 'GEN_PLOTS',1);

% for coneOS
params.UNDET_TOL = 1e-8;

% optvals taken from
% http://dimacs.rutgers.edu/Challenges/Seventh/Instances/tablestat.html
% http://people.orie.cornell.edu/miketodd/tttmpbc.pdf
objval = struct(...
    'nql30', -0.9460, ...
    'nql60', -0.935423, ...
    'nql180', 0, ... % n/a
    'qssp30', -6.496675, ...
    'qssp60', -6.562696, ...
    'qssp180', 0, ... % n/a
    'nb', -0.05070309, ...
    'nb_L1', -13.01227, ...
    'nb_L2', -1.628972, ...
    'nb_L2_bessel', -0.102571, ...
    'sched_50_50_orig', 26673.00, ...
    'sched_100_50_orig', 181889.9, ...
    'sched_100_100_orig', 717367, ...
    'sched_200_100_orig', 141360.4464, ...
    'sched_50_50_scaled',7.852038, ...
    'sched_100_50_scaled', 67.166281, ...   % their table is wrong or something
    'sched_100_100_scaled', 27.331457, ...
    'sched_200_100_scaled', 51.812471 ...
    );


pat = {'Total solve time is (?<total>[\d\.]+)', 'KKT matrix factorization took (?<fact>[\d\.]+)'};
pat_c = 'Time taken: (?<total>[\d\.]+)';

iter_pat = {'(?<iter>[\d]+) \|'};
iter_pat_c = {'(?<iter>[\d]+)\|'};

for i = 1:length(tests),
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
    cone = struct('f', 0, 'q', K.q', 'l', K.l, 'k_soc', length(K.q),'s', []);
    
    if cvx_on
        % CVX - sedumi:
        [m,n] = size(data.A);
        ns_soc = cone.q;
        idxs = [cone.f+cone.l; cone.f+cone.l + cumsum(ns_soc)];
        cvx_begin quiet
        cvx_solver sedumi
        variables xcvx(n) scvx(m)
        minimize (data.c'*xcvx)
        data.A*xcvx + scvx == data.b
        scvx(1:cone.f) == 0
        scvx(cone.f+1:cone.f+cone.l) >= 0
        for kk =1:cone.k_soc
            norm(scvx(idxs(kk) + 2: idxs(kk) + ns_soc(kk))) <= scvx(idxs(kk)+1)
        end
        cvx_end
        cvx_objval.(test_name) =  cvx_optval;
    end
    
    [output,x_m,z_m,status_m] = evalc('coneos_direct(data,cone,params);');
    
    coneos_err.(test_name) = 100 * abs( data.c'*x_m + objval.(test_name) ) / abs(objval.(test_name) + eps);
    
    timing = regexp(output, pat_c, 'names');
    coneos_times.(test_name) = str2num(timing.total);
    tmp = regexp(output,iter_pat_c, 'names');
    coneos_iters.(test_name) = str2num(tmp{1}(end).iter) + 1;
    coneos_objval.(test_name) = data.c'*x_m;
    if cvx_on
        coneos_perr.(test_name) = norm(x_m - xcvx);
        coneos_derr.(test_name) = norm(z_m - scvx);
    end
    [output,x_m,s_m,z_m,status_m] = evalc('pdos_direct(data,cone,params);');
    
    pdos_err.(test_name) = 100 * abs( data.c'*x_m + objval.(test_name) ) / abs(objval.(test_name) + eps);
    
    timing = regexp(output, pat, 'names');
    pdos_times.(test_name) = str2num(timing{1}.total);
    tmp = regexp(output,iter_pat, 'names');
    pdos_iters.(test_name) = str2num(tmp{1}(end).iter) + 1;
    pdos_objval.(test_name) = data.c'*x_m;
    if cvx_on
        pdos_perr.(test_name) = norm(x_m - xcvx);
        pdos_derr.(test_name) = norm(z_m - scvx);
    end
    
end

t_names = fieldnames(coneos_times);
for i = 1:numel(t_names)
    tdiff(i,1) = 100*(coneos_times.(t_names{i}) - pdos_times.(t_names{i}))/min(coneos_times.(t_names{i}), pdos_times.(t_names{i}));
    %objdiff(i,1) = 100*(coneos_objval.(t_names{i}) - pdos_objval.(t_names{i}))/min(coneos_objval.(t_names{i}),pdos_objval.(t_names{i}));
    errdiff(i,1) = coneos_err.(t_names{i}) - pdos_err.(t_names{i});
    if cvx_on
        pdiff(i,1) = 100*(coneos_perr.(t_names{i}) - pdos_perr.(t_names{i}))/min(coneos_perr.(t_names{i}),pdos_perr.(t_names{i}));
        ddiff(i,1) = 100*(coneos_derr.(t_names{i}) - pdos_derr.(t_names{i}))/min(coneos_derr.(t_names{i}),pdos_derr.(t_names{i}));
    end
end

save data/dimacs_run_data