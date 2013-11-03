#include "cones.h"
#ifdef LAPACK_LIB_FOUND
#include <cblas.h>
#include <lapacke.h>
#endif
void projectsdc(double * X, int n, Work * w); 

/* in place projection (with branches) */
void projCone(double *x, Cone * k, Work * w)
{
	int i;
	int count;

	/* project onto positive orthant */
	for(i = k->f; i < k->f+k->l; ++i)
	{
		if(x[i] < 0.0) x[i] = 0.0;
		//x[i] = (x[i] < 0.0) ? 0.0 : x[i];
	}
	count = k->l+k->f;
	/* project onto SOC */
	for(i = 0; i < k->qsize; ++i)
	{
		double v1 = x[count];
		double s = calcNorm(&(x[count+1]),k->q[i]-1);
		double alpha = (s + v1)/2.0;

		if(s <= v1) { /* do nothing */ }
		else if (s <= - v1) {
			memset(&(x[count]), 0, k->q[i]*sizeof(double));
		} else {    
			x[count] = alpha;
			scaleArray(&(x[count+1]), alpha/s, k->q[i]-1);
      //cblas_dscal(k->q[i]-1, alpha/s, &(x[count+1]),1);
		}           
		count += k->q[i];
	}
#ifdef LAPACK_LIB_FOUND
	/* project onto PSD cone */
	for (i=0; i < k->ssize; ++i){
		projectsdc(&(x[count]),k->s[i],w);
		count += (k->s[i])*(k->s[i]);
	}
#else
  if(k->ssize > 0){
    coneOS_printf("WARNING: solving SDP, no lapack library specified in makefile!\n");
    coneOS_printf("ConeOS will return a wrong answer!\n");
}
#endif
	/* project onto OTHER cones */
}

#ifdef LAPACK_LIB_FOUND
void projectsdc(double *X, int n, Work * w)
{ /* project onto the positive semi-definite cone */
  if (n == 1) {
    if(X[0] < 0.0) X[0] = 0.0; 
    return;
  }

  int i, j, m=0;
  double * Xs = w->Xs;
  double * Z = w->Z;
  double * e = w->e;
  memcpy(Xs,X,n*n*sizeof(double));

  // Xs = X + X', save div by 2 for eigen-recomp
  for (i = 0; i < n; ++i){
    cblas_daxpy(n, 1, &(X[i]), n, &(Xs[i*n]), 1);
    //b_daxpy(n, 1, &(X[i]), n, &(Xs[i*n]), 1);
  }
 
  double EIG_TOL = 1e-8;
  double vupper = calcNorm(Xs,n*n);
  LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'V', 'U', n, Xs, n, 0.0, vupper, -1, -1, EIG_TOL, &m, e, Z, n , NULL);
  //printf("m is %i, n is %i\n", m ,n);
  //printf("vupper is %f, max eig is %f\n",vupper, e[m>0 ? m-1:0]/2);
  memset(X, 0, n*n*sizeof(double));
  for (i = 0; i < m; ++i) {
    cblas_dsyr(CblasColMajor, CblasLower, n, e[i]/2, &(Z[i*n]), 1, X, n);
    //b_dsyr('L', n, -e[i]/2, &(Z[i*n]), 1, Xs, n);
  }
  // fill in upper half 
  for (i = 0; i < n; ++i){   
    for (j = i+1; j < n; ++j){   
      X[i + j*n] = X[j + i*n];    
    }   
  }
}
#endif

