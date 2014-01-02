function compileIndirect(flags)

common_coneOS = '../linAlg.c ../cones.c ../cs.c ../util.c ../coneOS.c coneOS_mex.c';
if (~isempty (strfind (computer, '64')))
    d = '-fPIC';
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end

% compile indirect
if (flags.COMPILE_WITH_OPENMP)
    cmd = sprintf('mex -O %s COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -std=c99 -O3 -fopenmp -pthread -DMATLAB_MEX_FILE %s %s" ../indirect/private.c %s -I../ %s -output coneos_indirect "LINKLIBS="\\$LINKLIBS -lm  %s %s"',  arr, d, flags.LCFLAG, common_coneOS, flags.INCS, flags.LOCS, flags.BLASLIB);
else
    cmd = sprintf('mex -O %s CFLAGS="-std=c99 -O3 -pthread -DMATLAB_MEX_FILE %s %s" ../indirect/private.c %s -I../ %s -output coneos_indirect LINKLIBS="\\$LINKLIBS -lm  %s %s"',  arr, d, flags.LCFLAG, common_coneOS, flags.INCS, flags.LOCS, flags.BLASLIB);
end
eval(cmd);
