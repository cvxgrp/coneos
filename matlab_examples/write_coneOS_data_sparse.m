function write_coneOS_data_sparse(data,K,params,name)

MAX_ITERS = 2500; % maximum num iterations for admm
EPS_ABS   = 1e-3; % quitting tolerances
UNDET_TOL = 1e-9; % tol for undetermined solution (tau = kappa = 0)
alpha=1.8;        % relaxation parameter (alpha = 1 is unrelaxed)
NORMALIZE = 1;
VERBOSE = 1;

% conjugate gradient (CG) settings:
USE_INDIRECT = false; % use conjugate gradient rather than direct method
CG_MAX_ITS = 15; % max iterations for CG
CG_TOL = 1e-9; % max CG quitting tolerance
CG_VERBOSE = false; % CG prints summary
%%
if ~isfield(params,'MAX_ITERS');params.MAX_ITERS = MAX_ITERS;end
if ~isfield(params,'EPS_ABS');params.EPS_ABS = EPS_ABS;end
if ~isfield(params,'UNDET_TOL');params.UNDET_TOL = UNDET_TOL;end
if ~isfield(params,'ALPHA');params.ALPHA = alpha;end
if ~isfield(params,'NORMALIZE');params.NORMALIZE = NORMALIZE;end
if ~isfield(params,'VERBOSE');params.VERBOSE = VERBOSE;end

% CG:
if ~isfield(params,'USE_INDIRECT');params.USE_INDIRECT = USE_INDIRECT;end
if ~isfield(params,'CG_MAX_ITS');params.CG_MAX_ITS = CG_MAX_ITS;end
if ~isfield(params,'CG_TOL');params.CG_TOL = CG_TOL;end
if ~isfield(params,'CG_VERBOSE');params.CG_VERBOSE = CG_VERBOSE;end

n = length(data.c);
m = size(data.A,1);
data.A = sparse(data.A);

%{
% symmetrize SD cone matrices:
[mm,nn]=size(data.A);
idx = K.f + K.l + sum(K.q);
for i=1:size(K.s)
    for j=1:nn
        work = data.A(idx+1:idx+K.s(i)^2, j);
        work = sparse(reshape(work,K.s(i),K.s(i)));
        if any(any(work~=work'))
            %fprintf('warning: symmetrizing A\n')
            work = (work+work')/2;
            data.A(idx+1:idx+K.s(i)^2, j) = sparse(work(:));
        end
    end
    
    work = data.b(idx+1:idx+K.s(i)^2);
    work = sparse(reshape(work,K.s(i),K.s(i)));
    if any(any(work~=work'))
        %fprintf('warning: symmetrizing b\n')
        work = (work+work')/2;
        data.b(idx+1:idx+K.s(i)^2) = sparse(work(:));
    end
    
    idx = idx + K.s(i)^2;
end
clear work;
%}
% col-compressed A
[i,~,s] = find(sparse(data.A));
i = i-1;
tmp = full(sum(data.A~=0));
pw = [0 cumsum(tmp)];
NNZ=length(i);

%Q=sparse([zeros(n) data.A' data.c;
%    -data.A zeros(m,m) data.b;
%    -data.c' -data.b' 0]);
%W=sparse([speye(n+m+1) Q';Q -speye(n+m+1)]);

delete(name);
fi = fopen(name,'w');
fprintf(fi,'%u ',n);fprintf(fi,'%u ',m);
fprintf(fi,'\n');
fprintf(fi,'%u ',K.f);fprintf(fi,'%u ',K.l);fprintf(fi,'%u ',length(K.q));fprintf(fi,'%u ',length(K.s));fprintf(fi,'%u ',K.ep);fprintf(fi,'%u ',K.ed);
fprintf(fi,'\n');
fprintf(fi,'%u ',params.MAX_ITERS);fprintf(fi,'%u ',params.CG_MAX_ITS);
fprintf(fi,'\n');
fprintf(fi,'%u ',params.VERBOSE); fprintf(fi,'%u ',params.NORMALIZE);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.ALPHA);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.UNDET_TOL);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.EPS_ABS);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.CG_TOL);
fprintf(fi,'\n');
fprintf(fi,'%u ',NNZ);
fprintf(fi,'\n');

fprintf(fi,'%u ',K.q');
fprintf(fi,'\n');
fprintf(fi,'%u ',K.s');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',full(data.b)');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',full(data.c)');
fprintf(fi,'\n');

fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',pw);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',s');
fprintf(fi,'\n');

% triplet A
%{
[i,j,s] = find(data.A);
i = i-1;j=j-1;
NNZ=length(i);
fprintf(fi,'%u ',NNZ);
fprintf(fi,'\n');
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',j');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',s');
fprintf(fi,'\n');
%}

fclose(fi);
