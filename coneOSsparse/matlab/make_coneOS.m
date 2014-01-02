% WARNING: OPENMP WITH MATLAB CAN CAUSE ERRORS AND CRASH, USE WITH CAUTION:
% openmp parallelizes the matrix multiply for the indirect solver (using CG):
flags.COMPILE_WITH_OPENMP = false;

% try to use blas libs; if not available then coneOS cannot solve SDPs:
try
    
    % EDIT THESE TO POINT TO YOUR BLAS + LAPACK LIBS:
    flags.BLASLIB = '-lopenblas -llapack -llapacke';
    flags.INCS = '-I/opt/local/include -I/usr/local/include';
    flags.LOCS = '-L/opt/local/lib -L/usr/local/lib -L/usr/lib';
    flags.LCFLAG = '-DLAPACK_LIB_FOUND';
    
    compileDirect(flags);
    compileIndirect(flags);
    
catch err
    
    flags.BLASLIB = '';
    flags.INCS = '';
    flags.LOCS = '';
    flags.LCFLAG = '';
    
    compileDirect(flags);
    compileIndirect(flags);    
    
    disp('Compiled without lapack support - unable to solve SDPs (can solve LPs, QPs, SOCPs, EXPs)')
    disp('To solve SDPs you must install cblas + lapacke and point the flags to the right locations')
end
