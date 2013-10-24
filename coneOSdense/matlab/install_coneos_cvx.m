clear all
make_coneOS

[ ver, isoctave, fs, ps ] = cvx_version;
coneos_install_path = strcat(fs,'/coneos');
mkdir(coneos_install_path);

%coneos_bin_loc = sprintf('%s/',coneos_install_path);
copyfile('coneos_direct.m*',coneos_install_path);
copyfile('coneos_indirect.m*',coneos_install_path);

copyfile('cvx_coneos.m', strcat(fs,'/shims/cvx_coneos.m'));
copyfile('coneos.m', coneos_install_path);
cvx_setup