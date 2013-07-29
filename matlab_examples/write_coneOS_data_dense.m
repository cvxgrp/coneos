function write_coneOS_data_dense(data,cone,params)

n = length(data.c);
m = size(data.A,1);
%Q=sparse([zeros(n) data.A' data.c;
%    -data.A zeros(m,m) data.b;
%    -data.c' -data.b' 0]);
%W=sparse([speye(n+m+1) Q';Q -speye(n+m+1)]);

delete data_dense;
fi = fopen('data_dense','w');
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


% dense A
fullA = full(data.A);
fprintf(fi,'%6.18f ',(fullA(:))');
fclose(fi);
