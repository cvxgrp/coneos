% the option -Dprintf=mexPrintf redfines all calls to "printf" with
%   "mexPrintf"
% switch(computer)
%     case {'PCWIN', 'GLNX86', 'MACI'}
%         evalc('system(''make -C .. purge packages CFLAGS="-m32 -DMATLAB_MEX_FILE"'')')
%     case {'PCWIN64', 'GLNXA64', 'SOL64', 'MACI64'}
%         evalc('system(''make -C .. purge packages CFLAGS="-DMATLAB_MEX_FILE"'')')
%     otherwise
%         evalc('system(''make -C .. purge packages CFLAGS="-DMATLAB_MEX_FILE"'')')
% end

%compile direct
common_coneOS = 'coneOS_mex.c ../coneOS.c ../cones.c ../cs.c ../linAlg.c ../util.c';


% link your blas libraries
% gsl blas:
%BLASLIB = '-lgslcblas';

% open blas:
BLASLIB = '-lopenblas';

% multi-threaded atlas:
%BLASLIB = '-L/usr/local/atlas/lib -lpthread -lptcblas -lptf77blas -llapack -latlas';

% single-threaded atlas:
%BLASLIB = '-L/usr/local/atlas/lib -lcblas -lf77blas -llapack -latlas';

if (~isempty (strfind (computer, '64')))
    d = '-fPIC' ;
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end

cmd = sprintf ('mex -v -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s" -I../', arr, d) ;
amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
    'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
for i = 1 : length (amd_files)
    cmd = sprintf ('%s ../direct/%s.c', cmd, amd_files {i}) ;
end
cmd = sprintf ('%s ../direct/ldl.c %s ../direct/private.c -lm -o coneOS_direct', cmd, common_coneOS) ;
eval(cmd) ;

cmd = sprintf ('mex -v -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s" -I../', arr, d) ;
% compile GSL direct
cmd = sprintf ('%s ../direct/ldl.c %s ../direct_gsl/private.c -I/opt/local/include -L/opt/local/lib -lm -lgsl %s -o coneOS_direct_gsl', cmd, common_coneOS,BLASLIB) ;
eval(cmd) ;

% compile indirect (XXX: openmp?)
cmd = sprintf('mex -v -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s" %s ../indirect/private.c -I../ -o coneOS_indirect -lm', arr, d, common_coneOS);
eval(cmd);

% compile GSL indirect
cmd = sprintf('mex -v -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s" %s ../indirect_gsl/private.c -I../ -I/opt/local/include -L/opt/local/lib -o coneOS_indirect_gsl -lm -lgsl %s', arr, d, common_coneOS,BLASLIB);
eval(cmd) ;
