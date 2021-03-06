#include "mex.h" 
#include "matrix.h"
#include "coneOS.h"
#include "linAlg.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* matlab usage: coneOS(data,cone,params); */
  
  if (nrhs != 3){
    mexErrMsgTxt("Three arguments are required in this order: data struct, cone struct, params struct");
  }
  Data * d = mxMalloc(sizeof(Data)); 
  Cone * k = mxMalloc(sizeof(Cone)); 
  const mxArray *data = prhs[0];
   
  const mxArray *A_mex = (mxArray *) mxGetField(data,0,"A");
  if(A_mex == NULL) {
    coneOS_free(d); coneOS_free(k);
    mexErrMsgTxt("Data struct must contain a `A` entry.");
  }
  if (!mxIsSparse(A_mex)){
    coneOS_free(d); coneOS_free(k);
    mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
  }
  const mxArray *b_mex = (mxArray *) mxGetField(data,0,"b");
  if(b_mex == NULL) {
    coneOS_free(d); coneOS_free(k);
    mexErrMsgTxt("Data struct must contain a `b` entry.");
  }
  const mxArray *c_mex = (mxArray *) mxGetField(data,0,"c"); 
  if(c_mex == NULL) {
    coneOS_free(d); coneOS_free(k);
    mexErrMsgTxt("Data struct must contain a `c` entry.");
  }

  const mxArray *cone = prhs[1];
  const mxArray *params = prhs[2];
  d->n = *(mxGetDimensions(c_mex));
  d->m = *(mxGetDimensions(b_mex));
 
  d->b = mxGetPr(b_mex);
  d->c = mxGetPr(c_mex);
  
  const mxArray *ALPH_mex = mxGetField(params,0,"ALPHA");
  if (ALPH_mex == NULL) d->ALPH = 1.8;
  else d->ALPH = (double)*mxGetPr(ALPH_mex);
  
  const mxArray *RHO_X_mex = mxGetField(params,0,"RHO_X");
  if (RHO_X_mex == NULL) d->RHO_X = 1e-3;
  else d->RHO_X = (double)*mxGetPr(RHO_X_mex);

  const mxArray *MAX_ITERS_mex = mxGetField(params,0,"MAX_ITERS");
  if (MAX_ITERS_mex == NULL) d->MAX_ITERS = 2500;
  else d->MAX_ITERS = (int)*mxGetPr(MAX_ITERS_mex);

  const mxArray *EPS_ABS_mex = mxGetField(params,0,"EPS_ABS");
  if (EPS_ABS_mex == NULL) d->EPS_ABS = 1e-3;
  else d->EPS_ABS = (double)*mxGetPr(EPS_ABS_mex);

  const mxArray *UNDET_TOL_mex = mxGetField(params,0,"UNDET_TOL");
  if (UNDET_TOL_mex == NULL) d->UNDET_TOL = 1e-9;
  else d->UNDET_TOL = (double)*mxGetPr(UNDET_TOL_mex);

  const mxArray *VERBOSE_mex = mxGetField(params,0,"VERBOSE");
  if (VERBOSE_mex == NULL) d->VERBOSE = 0;
  else d->VERBOSE = (int)*mxGetPr(VERBOSE_mex);

  const mxArray *NORMALIZE_mex = mxGetField(params,0,"NORMALIZE");
  if (NORMALIZE_mex == NULL) d->NORMALIZE = 1;
  else d->NORMALIZE = (int)*mxGetPr(NORMALIZE_mex);

  k->f = (int)*mxGetPr(mxGetField(cone,0,"f"));
  k->l = (int)*mxGetPr(mxGetField(cone,0,"l"));
  
  mxArray * ep = mxGetField(cone,0,"ep");
  if(ep) k->ep = (int)*mxGetPr(ep);
  else k->ep = 0;
  mxArray * ed = mxGetField(cone,0,"ed");
  if(ed) k->ed = (int)*mxGetPr(ed);
  else k->ed = 0;
  
  double * q_mex = mxGetPr(mxGetField(cone,0,"q"));
  k->qsize = *(mxGetDimensions(mxGetField(cone,0,"q")));
  int i;
  k->q = mxMalloc(sizeof(int)*k->qsize);
  for ( i=0; i < k->qsize; i++ ){
    k->q[i] = (int)q_mex[i]; 
  }
  
  double * s_mex = mxGetPr(mxGetField(cone,0,"s"));
  k->ssize = *(mxGetDimensions(mxGetField(cone,0,"s")));
  k->s = mxMalloc(sizeof(int)*k->ssize);
  for ( i=0; i < k->ssize; i++ ){
    k->s[i] = (int)s_mex[i]; 
  }
  
  d->Anz = mxGetNzmax(A_mex);
  d->Ax = (double *)mxGetPr(A_mex);
  d->Ap = (int *)mxMalloc(sizeof(int)*d->Anz);
  d->Ai = (int *)mxMalloc(sizeof(int)*d->Anz);
  long * A_i = (long*)mxGetIr(A_mex);
  /* XXX fix this crap: */
  for (i = 0; i < d->Anz; i++){
    d->Ai[i] = (int)A_i[i];
  }
  long * A_p = (long*)mxGetJc(A_mex);
  for (i = 0; i < (d->n)+1; i++) {          
    d->Ap[i] = (int)A_p[i];
  }
  
  Sol sol;
  Info info;

  int status = coneOS(d,k,&sol,&info);

  plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
  mxSetPr(plhs[0], sol.x);
  mxSetM(plhs[0], d->n); 
  mxSetN(plhs[0], 1); 

  plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
  mxSetPr(plhs[1], sol.y);
  mxSetM(plhs[1], d->m); 
  mxSetN(plhs[1], 1); 

  const char * infoFields[] = {"iter","status","pobj","dobj","resPri","resDual","relGap","time"}; 
  const int numInfoFields = 8;
  mwSize one[1] = {1};
  mxArray * xm;
  plhs[2] = mxCreateStructArray(1,one,numInfoFields,infoFields);

  mxSetField(plhs[2], 0, "status", mxCreateString(info.status));
   
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[2], 0, "iter", xm);
  *mxGetPr(xm) = info.iter;
  
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[2], 0, "pobj", xm);
  *mxGetPr(xm) = info.pobj;

  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[2], 0, "dobj", xm);
  *mxGetPr(xm) = info.dobj;
  
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[2], 0, "resPri", xm);
  *mxGetPr(xm) = info.resPri;
  
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[2], 0, "resDual", xm);
  *mxGetPr(xm) = info.resDual;
  
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[2], 0, "relGap", xm);
  *mxGetPr(xm) = info.relGap;

  //info.time is millisecs - return value in secs
  xm = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxSetField(plhs[2], 0, "time", xm);
  *mxGetPr(xm) = info.time/1e3; 

  if(k->q) coneOS_free(k->q);
  if(k->s) coneOS_free(k->s);
  if(d->Ap) coneOS_free(d->Ap);
  if(d->Ai) coneOS_free(d->Ai);
  if(d) coneOS_free(d);
  if(k)  coneOS_free(k);
  return; 
}

