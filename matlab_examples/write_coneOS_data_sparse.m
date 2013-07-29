function write_coneOS_data_sparse(data,cone,params)

n = length(data.c);
m = size(data.A,1);

%Q=sparse([zeros(n) data.A' data.c;
%    -data.A zeros(m,m) data.b;
%    -data.c' -data.b' 0]);
%W=sparse([speye(n+m+1) Q';Q -speye(n+m+1)]);

delete data_sparse;
fi = fopen('data_sparse','w');
fprintf(fi,'%u ',n);fprintf(fi,'%u ',m);
fprintf(fi,'\n');
fprintf(fi,'%u ',cone.f);fprintf(fi,'%u ',cone.l);fprintf(fi,'%u ',length(cone.q));fprintf(fi,'%u ',length(cone.s));
fprintf(fi,'\n');
fprintf(fi,'%u ',cone.q');
fprintf(fi,'\n');
fprintf(fi,'%u ',cone.s');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',data.b');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',data.c');
fprintf(fi,'\n');

fprintf(fi,'%u ', params.MAX_ITERS);fprintf(fi,'%u ',params.CG_MAX_ITS);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.ALPHA);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.UNDET_TOL);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.EPS_ABS);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',params.CG_TOL);
fprintf(fi,'\n');
fprintf(fi,'%u ',params.VERBOSE); fprintf(fi,'%u ',params.NORMALIZE);
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

% col-compressed A
[i,j,s] = find(sparse(data.A));
i = i-1;
tmp = full(sum(data.A~=0));
pw = [0 cumsum(tmp)];
NNZ=length(i);
fprintf(fi,'%u ',NNZ);
fprintf(fi,'\n');
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',pw);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',s');
fprintf(fi,'\n');

%{
[i,j,s] = find(W);
i = i-1;
tmp = full(sum(W)~=0);
pw = [0 cumsum(tmp)];
NNZ=length(i);
fprintf(fi,'%u ',NNZ);
fprintf(fi,'\n');
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',pw);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',s');
fprintf(fi,'\n');
%}
fclose(fi);
