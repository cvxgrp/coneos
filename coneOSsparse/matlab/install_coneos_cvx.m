make_coneOS

cvx_path = which('cvx_setup.m')
cvx_path = cvx_path(1:end-length('cvx_setup.m'))
coneos_install_path = sprintf('%sconeos',cvx_path);
mkdir(coneos_install_path);
coneos_bin_loc = sprintf('%s/coneos.mexmaci64',coneos_install_path);
copyfile('coneOS_direct.mexmaci64',coneos_bin_loc);
copyfile('cvx_coneos.m', sprintf('%sshims/cvx_coneos.m',cvx_path));

cvx_setup