clear all

cvx_path = which('cvx_setup.m')
cvx_path = cvx_path(1:end-length('cvx_setup.m'))
coneos_install_path = sprintf('%sconeos_matlab',cvx_path);
mkdir(coneos_install_path);

%coneos_bin_loc = sprintf('%s/',coneos_install_path);
copyfile('coneos_matlab.m',coneos_install_path);

copyfile('cvx_coneos_matlab.m', sprintf('%sshims/cvx_coneos_matlab.m',cvx_path));
cvx_setup
