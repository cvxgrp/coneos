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

common_coneOS = 'coneOS_mex.c ../coneOS.c ../cones.c ../cs.c ../linAlg.c ../util.c';
if (~isempty (strfind (computer, '64')))
    d = '-fPIC';
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end

% try to use blas libs, if not available can only solve SOCPs:
try
   
    % EDIT THESE TO POINT TO YOUR BLAS + LAPACK LIBS:
    BLASLIB = '-lopenblas -llapack -llapacke';
    INCS = '-I/opt/local/include -I/usr/local/include';
    LOCS = '-L/opt/local/lib -L/usr/local/lib';
    LCFLAG = '-DLAPACK_LIB_FOUND';
        
    %compile direct
    cmd = sprintf ('mex -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s %s" -I../ %s', arr, d, LCFLAG, INCS) ;
    amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
        'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
        'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
    for i = 1 : length (amd_files)
        cmd = sprintf ('%s ../direct/%s.c', cmd, amd_files {i}) ;
    end
    cmd = sprintf ('%s ../direct/ldl.c %s ../direct/private.c -lm %s %s -o coneos_direct', cmd, common_coneOS, LOCS, BLASLIB) ;
    eval(cmd);
    
    % compile indirect (XXX: openmp?)
    %cmd = sprintf('mex -O %s COMPFLAGS="/openmp $COMPFLAGS" CFLAGS="\\$CFLAGS -std=c99 -O3 -fopenmp -pthread -DMATLAB_MEX_FILE %s %s" ../indirect/private.c %s -I../ %s -o coneos_indirect LDFLAGS="\\$LDFLAGS -fopenmp -lm  %s %s"',  arr, d, LCFLAG, common_coneOS, INCS, LOCS, BLASLIB);
    cmd = sprintf('mex -O %s CFLAGS="-std=c99 -O3 -pthread -DMATLAB_MEX_FILE %s %s" ../indirect/private.c %s -I../ %s -o coneos_indirect LDFLAGS="\\$LDFLAGS -lm  %s %s"',  arr, d, LCFLAG, common_coneOS, INCS, LOCS, BLASLIB);
    eval(cmd);
    %mex -v -O COMPFLAGS="/openmp $COMPFLAGS" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims ../coneOS.c ../linAlg.c ../cones.c ../cs.c ../util.c coneOS_mex.c ../indirect/private.c -I../ -output coneOS_indirect -lm
    
catch err
    
    disp('Compiling without lapack support - can only solve SOCPs')
    
    BLASLIB = '';
    INCS = '';
    LOCS = '';
    LCFLAG = '';
    
    %compile direct
    cmd = sprintf ('mex -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s %s" -I../ %s', arr, d, LCFLAG, INCS) ;
    amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
        'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
        'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
    for i = 1 : length (amd_files)
        cmd = sprintf ('%s ../direct/%s.c', cmd, amd_files {i}) ;
    end
    cmd = sprintf ('%s ../direct/ldl.c %s ../direct/private.c -lm %s %s -o coneos_direct', cmd, common_coneOS, LOCS, BLASLIB) ;
    eval(cmd);
    
    % compile indirect (XXX: openmp?)
    %cmd = sprintf('mex -O %s COMPFLAGS="/openmp $COMPFLAGS" CFLAGS="\\$CFLAGS -std=c99 -O3 -fopenmp -pthread -DMATLAB_MEX_FILE %s %s" ../indirect/private.c %s -I../ %s -o coneos_indirect LDFLAGS="\\$LDFLAGS -fopenmp -lm  %s %s"',  arr, d, LCFLAG, common_coneOS, INCS, LOCS, BLASLIB);
    cmd = sprintf('mex -O %s CFLAGS="-std=c99 -O3 -pthread -DMATLAB_MEX_FILE %s %s" ../indirect/private.c %s -I../ %s -o coneos_indirect LDFLAGS="\\$LDFLAGS -lm  %s %s"',  arr, d, LCFLAG, common_coneOS, INCS, LOCS, BLASLIB);
    eval(cmd);
    %mex -v -O COMPFLAGS="/openmp $COMPFLAGS" CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims ../coneOS.c ../linAlg.c ../cones.c ../cs.c ../util.c coneOS_mex.c ../indirect/private.c -I../ -output coneOS_indirect -lm
    
    
end


