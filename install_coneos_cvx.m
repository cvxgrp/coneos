path = which('cvx_setup.m')
path = path(1:end-length('cvx_setup.m'))
mkdir(sprintf('%sconeos',path));
