common_coneOS = 'coneOS_mex.c ../coneOS.c ../linAlg.c ../cones.c ../util.c';

% blas + lapack libraries:
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
cmd = sprintf ('mex -O %s CFLAGS="-std=c99 -pedantic -O3 -DMATLAB_MEX_FILE -DLAPACK_LIB_FOUND %s" -I../', arr, d) ;
cmd = sprintf ('%s %s ../direct/private.c -I../external/ %s -lm %s -output coneos_direct', cmd, common_coneOS, LOCS, BLASLIB) ;
eval(cmd);


% compile indirect
cmd = sprintf ('mex -O %s CFLAGS="-std=c99 -pedantic -O3 -DMATLAB_MEX_FILE -DLAPACK_LIB_FOUND %s" -I../', arr, d) ;
cmd = sprintf ('%s %s ../indirect/private.c -I../external/ %s -lm %s -output coneos_indirect', cmd, common_coneOS, LOCS, BLASLIB) ;
eval(cmd);
