common_coneOS = 'coneOS_mex.c ../coneOS.c ../cones.c ../linAlg.c ../util.c ../external/lapacke.c';


% open blas:
BLASLIB = '-lmwlapack -lcblas -lmwblas';
%BLASLIB = '-lopenblas';

if (~isempty (strfind (computer, '64')))
    d = '-fPIC' ;
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end
LOCS = '-I/opt/local/include -I/usr/local/include -L/opt/local/lib -L/usr/local/lib'

% compile GSL direct
cmd = sprintf ('mex -v -O %s CFLAGS="-std=c99 -pedantic -O3 -DMATLAB_MEX_FILE %s" -I../', arr, d) ;
cmd = sprintf ('%s %s ../direct/private.c -I../external/ %s -lm %s -o coneOS', cmd, common_coneOS, LOCS, BLASLIB) ;
eval(cmd) ;

%{
% compile indirect (XXX: openmp?)
cmd = sprintf('mex -v -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s" %s ../indirect/private.c -I../ -o coneOS_indirect -lm', arr, d, common_coneOS);
eval(cmd);

% compile GSL indirect
cmd = sprintf('mex -v -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s" %s ../indirect_gsl/private.c -I../ -I/opt/local/include -L/opt/local/lib -o coneOS_indirect_gsl -lm -lgsl %s', arr, d, common_coneOS,BLASLIB);
eval(cmd) ;
%}