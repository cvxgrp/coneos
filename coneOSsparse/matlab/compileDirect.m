function compileDirect(flags)

common_coneOS = '../linAlg.c ../cones.c ../cs.c ../util.c ../coneOS.c coneOS_mex.c';
if (~isempty (strfind (computer, '64')))
    d = '-fPIC';
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end

cmd = sprintf ('mex -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s %s" -I"../" %s', arr, d, flags.LCFLAG, flags.INCS) ;
amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
    'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
for i = 1 : length (amd_files)
    cmd = sprintf ('%s ../direct/%s.c', cmd, amd_files {i}) ;
end
cmd = sprintf ('%s ../direct/ldl.c %s ../direct/private.c LINKLIBS="\\$LINKLIBS -lm %s %s" -output coneos_direct', cmd, common_coneOS, flags.LOCS, flags.BLASLIB) ;
eval(cmd);
