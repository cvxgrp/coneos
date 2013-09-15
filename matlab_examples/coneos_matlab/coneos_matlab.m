function [x,z,info] = coneos_matlab(data,K,params)
% cone solver, solves:
%
% min. c'x
% subject to Ax + s = b
% s \in K
%
% where x \in R^n, s \in R^m
%
% K is product of cones in this particular order:
% zero cone, lp cone, second order cone(s)
%
% data must consist of data.A, data.b, data.c, where A,b,c used as above
%
% cone struct must consist of:
% cone.zero, length of zero cone (for equality constraints)
% cone.lp, length of lp cone
% cone.k_soc, number of second-order cones (SOC)
% cone.q, array of SOC lengths
% (thus m = cones.zero+cone.lp+sum(cone.q) = size(A,1))
%
% cone struct is only used in proj_cone, to add other cones
% simply add the relevant size data to the cone struct and edit the
% proj_cone method to include projection onto the new cone AND
% the dual cone (current implementation is second-order cones only)

% params struct consists of the following fields
% (here set to default settings):
GEN_PLOTS = false;% generate convergence plots
MAX_ITERS = 2000; % maximum num iterations for admm
EPS_ABS   = 5e-4; % quitting tolerances
UNDET_TOL = 1e-9; % tol for undetermined solution (tau = kappa = 0)
alpha=1.8;        % relaxation parameter (alpha = 1 is unrelaxed)
NORMALIZE = 1;
RELAX_X = false;

% conjugate gradient (CG) settings:
USE_INDIRECT = false; % use conjugate gradient rather than direct method
CG_MAX_ITS = 30; % max iterations for CG
CG_TOL = 1e-9; % max CG quitting tolerance
CG_VERBOSE = false; % CG prints summary
% experimental:
rho_x = 1e-3;
%%
if nargin==3
    if isfield(params,'GEN_PLOTS');GEN_PLOTS = params.GEN_PLOTS;end
    if isfield(params,'MAX_ITERS');MAX_ITERS = params.MAX_ITERS;end
    if isfield(params,'EPS_ABS');EPS_ABS = params.EPS_ABS;end
    if isfield(params,'UNDET_TOL');UNDET_TOL = params.UNDET_TOL;end
    if isfield(params,'ALPHA');alpha = params.ALPHA;end
    if isfield(params,'NORMALIZE');NORMALIZE = params.NORMALIZE;end
    % experimental
    if isfield(params,'RHO_X');rho_x = params.RHO_X;end
    if isfield(params,'RELAX_X');RELAX_X = params.RELAX_X;end
    % CG:
    if isfield(params,'USE_INDIRECT');USE_INDIRECT = params.USE_INDIRECT;end
    if isfield(params,'CG_MAX_ITS');CG_MAX_ITS = params.CG_MAX_ITS;end
    if isfield(params,'CG_TOL');CG_TOL = params.CG_TOL;end
    if isfield(params,'CG_VERBOSE');CG_VERBOSE = params.CG_VERBOSE;end
end

%%

n = length(data.c);
m = length(data.b);
l=n+m+1;

%{
Q=sparse([zeros(n) data.A' data.c;
    -data.A zeros(m,m) data.b;
   -data.c' -data.b' 0]);
%}
if (CG_MAX_ITS > l)
    CG_MAX_ITS = l;
end

%%

if (NORMALIZE)
    disp('ORIGINAL NORMALIZATION SCHEME')
    A_o = data.A;
    b_o = data.b;
    c_o = data.c;
    D = ones(m,1);
    E = ones(n,1);
    NN = 1; % NN = 1, other choices bad
    for j=1:NN
        % D scale:
        Dt = norms(data.A(1:K.f,:)')';
        idx = K.f;
        Dt = [Dt;norms(data.A(idx+1:idx+K.l,:)')'];
        idx = idx + K.l;
        for i=1:length(K.q)
            nmA = mean(norms(data.A(idx+1:idx+K.q(i),:)'));
            Dt = [Dt;nmA*ones(K.q(i),1)];
            idx = idx + K.q(i);
        end
        for i=1:length(K.s)
            nmA = mean(norms(data.A(idx+1:idx+K.s(i)^2,:)'));
            Dt = [Dt;nmA*ones(K.s(i)^2,1)];
            idx = idx + K.s(i)^2;
        end
        D = max(D.*Dt,1e-6);
        data.A = sparse(diag(1./Dt))*data.A;
        % E Scale
        Et = max(norms(data.A),1e-6)';
        E = E.*Et;
        data.A = data.A*sparse(diag(1./Et));
    end
    nmrowA = mean(norms(data.A'));
    
    data.b = data.b./D;
    sc_b = 1/max(norm(data.b),1e-4);
    data.b = full(data.b*sc_b);
    
    data.c = data.c./E;
    sc_c = nmrowA/max(norm(data.c),1e-4);
    data.c = data.c*sc_c;

    rr = 4;
    data.A=rr*data.A;
    data.c=rr*data.c;
    data.b=rr*data.b;
    %{
    norm(D)
    norm(E)
    norm(norms(data.A))
    norm(data.b)
    norm(data.c)
    %}
end

%%
work = struct('USE_INDIRECT', USE_INDIRECT);
if ~USE_INDIRECT
    W=sparse([rho_x*speye(n) data.A';data.A -speye(m)]);
    disp('Factorization')
    try
        work.P=amd(W);
        [work.L,work.D] = ldlsparse(W,work.P);
        work.L=work.L+speye(n+m);
    catch ldlerror
        disp('WARNING: LDLSPARSE ERROR, using MATLAB LDL instead (this is slower).')
        [work.L,work.D,work.P] = ldl(W,'vector');
    end
else
    work.CG_VERBOSE = CG_VERBOSE;
end

h = [data.c;data.b];
[g,itn] = solveLinSystem(work,data,h,n,m,CG_MAX_ITS*100,CG_TOL/100,zeros(n+m,1),rho_x);
g(n+1:end) = -g(n+1:end);
gTh = g'*h;


% u = [x;z;tau], v = [y;s;kappa]
disp('Solving cone program')

if (isfield(params,'WARMXY'))
    u = [params.WARMXY;1];
    v = [zeros(n,1);data.b*u(end) - data.A*u(1:n);0];
    ut = u;
else
    u = zeros(l,1);u(end) = 1;
    v = zeros(l,1);%v(end) = 0;
end

%G = eye(n) + data.A'*data.A;

for i=1:MAX_ITERS
    u_prev = u;
    % solve linear system
    ut = u+v;
    ut(1:n) = rho_x*ut(1:n);
    ut(1:n+m) = ut(1:n+m) - ut(end)*h;
    ut(1:n+m) = ut(1:n+m) - h*((g'*ut(1:n+m))/(gTh+1));
    warm_start = u(1:n+m);
    ut(n+1:end-1) = -ut(n+1:end-1);
    %%
    [ut(1:n+m),itn] = solveLinSystem(work,data,ut(1:n+m),n,m,CG_MAX_ITS,CG_TOL,warm_start,rho_x);
    ut(end) = (ut(end) + h'*ut(1:n+m));
    
    % K proj:
    rel_ut = alpha*ut+(1-alpha)*u;
    if (~RELAX_X)
        rel_ut(1:n) = ut(1:n);
    end
    u = rel_ut - v;
    u(n+1:n+m) = proj_cone(u(n+1:n+m),K);
    u(l) = pos(u(l));
    
    % dual update:
    v = v + (u - rel_ut);
    %vt = v + u - u_old;
    
    %% convergence checking:
    tau = 0.5*abs(u(l)+ut(l));%ut(l);
    kap = abs(v(end));
    
    nm = 2;
    err_pri = norm(u-ut,nm);%/(tau+kap);
    err_dual = norm(u-u_prev,nm);%/(tau+kap);
    if GEN_PLOTS
        if USE_INDIRECT; cg_its(i) = itn; mults(i) = 2+2*itn;end
        nms(i,1) = err_pri;
        nms(i,2) = err_dual;
        tau_i(i) = ut(end);
        kap_i(i) = v(end);
        pobj(i) = data.c'*ut(1:n)/tau;
        dobj(i) = -data.b'*ut(n+1:n+m)/tau;
        %utv(i) = abs(ut'*v/(tau+kap));
        %vtu(i) = abs(vt'*u/(tau+kap));
    end
    
    if (min(tau,kap)/max(tau,kap) < 1e-6 && max(err_pri,err_dual) < EPS_ABS*(tau+kap))
        break
    end
    
    if mod(i-1,100)==0
        fprintf('Iteration %i, primal residual %4e, dual residual %4e, kap/tau %4e\n',i-1,err_pri/(tau+kap),err_dual/(tau+kap),kap/tau);
    end
end
fprintf('Iteration %i, primal residual %4e, dual residual %4e, kap/tau %4e\n',i-1,err_pri/(tau+kap),err_dual/(tau+kap),kap/tau);

%%
tau = 0.5*(u(l)+ut(l));
kap = v(end);
x_h = u(1:n);
z_h = 0.5*(u(n+1:n+m) + ut(n+1:n+m));
%s_h = v(n+1:n+m);
% tau, kap checking not exactly consistent w/ Lieven notes
% since we allow small numerical errors in solution
if (tau > UNDET_TOL && tau > kap) % this is different to Lieven
    status = 'Solved'
    x=x_h/tau;
    z=z_h/tau;
    %s=s_h/tau;
else
    x = nan(n,1);
    z = nan(m,1);
    if norm((u+ut)/2)<=2*(UNDET_TOL*sqrt(l))
        status = 'Undetermined'
    elseif data.b'*z_h < data.c'*x_h % this is different to Lieven
        status = 'Infeasible'
        z = -z_h/(data.b'*z_h);
    else
        status = 'Unbounded'
        x = -x_h/(data.c'*x_h);
    end
end
fprintf('Took %i iterations\n',i)
info.status = status;
info.iter = i;

if (NORMALIZE)
    z = z./(D*sc_c);
    x = x./(E*sc_b);
end

%%
if GEN_PLOTS
    figure();semilogy(nms(:,1));hold on;semilogy(nms(:,2),'r');
    legend('pri resid','dual resid');
    figure();plot(tau_i);hold on; plot(kap_i,'r');
    legend('tau','kappa')
    figure();plot(pobj);hold on;plot(dobj,'r');
    legend('primal obj','dual obj')
    %figure(); semilogy(utv);hold on;semilogy(vtu,'r');
    if USE_INDIRECT;
        figure();plot(cg_its);xlabel('k');ylabel('Conjugate Gradient Iterations');
        figure();semilogy(cumsum(mults),nms(:,1));hold on;semilogy(cumsum(mults),nms(:,2),'r');
        xlabel('A multiplies');ylabel('norm error');
    end
end
end


function z = proj_cone(z,c)
free_len = c.f;
lp_len = c.l;
k_soc = length(c.q);
q = c.q;
s = c.s;
ssize = length(c.s);
% lp cone
z(free_len+1:lp_len+free_len) = pos(z(free_len+1:lp_len+free_len));
% SOCs
idx=lp_len+free_len;
for i=1:k_soc
    z(idx+1:idx+q(i)) = proj_soc(z(idx+1:idx+q(i)));
    idx=idx+q(i);
end
%SDCs
for i=1:ssize
    z(idx+1:idx+s(i)^2) = proj_sdp(z(idx+1:idx+s(i)^2),s(i));
    idx=idx+s(i)^2;
end
end

function z = proj_soc(tt)
v1=tt(1);v2=tt(2:end);
if norm(v2)<=-v1
    v2=zeros(length(v2),1);v1=0;
elseif norm(v2)> abs(v1)
    v2=0.5*(1+v1/norm(v2))*v2;
    v1=norm(v2);
end
z=[v1;v2];
end

function z = proj_sdp(z,n)
z = reshape(z,n,n);
zs=z+z';
[V,S] = eig(zs);
S(S>0) = 0;
z = V*(-S/2)*V' + z;
z = z(:);
end

function [y,itn] = solveLinSystem(w,data,rhs,n,m,CG_MAX_ITERS,CG_TOL,warm_start,rho_x)
% assumes matrix [I A';A -I]
if w.USE_INDIRECT
    y=zeros(n+m,1);
    [y(1:n), itn] = pcg_A(data.A,rhs(1:n)+data.A'*rhs(n+1:n+m),warm_start(1:n),rho_x,CG_MAX_ITERS,CG_TOL,w.CG_VERBOSE);
    y(n+1:n+m) = -rhs(n+1:n+m) + data.A*y(1:n);
else
    % SOLVE DIRECTLY:
    y = (w.L'\(w.D\(w.L\rhs(w.P))));y(w.P)=y;
    itn = -1;
end
end

function [x,i] = pcg_A(A,b,x,rho_x,MAX_ITS,TOL,VERBOSE)
r=b-(rho_x*x+A'*(A*x));
p=r;
rsold=r'*r;

for i=1:MAX_ITS
    Ap=rho_x*p + A'*(A*p);
    alpha=rsold/(p'*Ap);
    x=x+alpha*p;
    r=r-alpha*Ap;
    rsnew=r'*r;
    resid = sqrt(rsnew);
    if resid<TOL
        if VERBOSE
            fprintf('CG took %i iterations to converge, resisdual %4f <= tolerance %4f\n',i,resid,TOL)
        end
        return;
    end
    p=r+(rsnew/rsold)*p;
    rsold=rsnew;
end
if VERBOSE
    fprintf('CG did not converge within %i iterations, resisdual %4f > tolerance %4f\n',MAX_ITS,resid,TOL)
end
end