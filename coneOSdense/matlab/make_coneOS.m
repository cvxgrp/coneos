common_coneOS = 'coneOS_mex.c ../coneOS.c ../linAlg.c ../cones.c ../util.c';
cd ..
evalc('system(''make purge'')');
cd matlab

% open blas:
%BLASLIB = '-lmwlapack -lcblas -lmwblas';
BLASLIB = '-lopenblas -llapack -llapacke';

if (~isempty (strfind (computer, '64')))
    d = '-fPIC' ;
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end
LOCS = '-I/opt/local/include -I/usr/local/include -L/opt/local/lib -L/usr/local/lib';

% compile direct
delete('coneOS_direct.mex*')
cmd = sprintf ('mex -O %s CFLAGS="-std=c99 -pedantic -O3 -DMATLAB_MEX_FILE -DLAPACK_LIB_FOUND %s" -I../', arr, d) ;
cmd = sprintf ('%s %s ../direct/private.c -I../external/ %s -lm %s -o coneos_direct', cmd, common_coneOS, LOCS, BLASLIB) ;
eval(cmd) ;


% compile indirect
cmd = sprintf ('mex -O %s CFLAGS="-std=c99 -pedantic -O3 -DMATLAB_MEX_FILE -DLAPACK_LIB_FOUND %s" -I../', arr, d) ;
cmd = sprintf ('%s %s ../indirect/private.c -I../external/ %s -lm %s -o coneos_indirect', cmd, common_coneOS, LOCS, BLASLIB) ;
eval(cmd) ;
