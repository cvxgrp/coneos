clear all
make_coneOS

cvx_path = which('cvx_setup.m')
cvx_path = cvx_path(1:end-length('cvx_setup.m'))
coneos_install_path = sprintf('%sconeos',cvx_path);
mkdir(coneos_install_path);

%coneos_bin_loc = sprintf('%s/',coneos_install_path);
copyfile('coneos_direct.mex*',coneos_install_path);
copyfile('coneos_indirect.mex*',coneos_install_path);

copyfile('cvx_coneos.m', sprintf('%sshims/cvx_coneos.m',cvx_path));
copyfile('coneos.m', coneos_install_path);
cvx_setup
