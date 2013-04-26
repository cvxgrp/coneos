#include "cones.h"
#include <cblas.h>
#include <lapacke.h>

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
		double s = cblas_dnrm2(k->q[i]-1, &(x[count+1]), 1);
		double alpha = (s + v1)/2.0;

		if(s <= v1) { /* do nothing */ }
		else if (s <= - v1) {
			memset(&(x[count]), 0, k->q[i]*sizeof(double));
		} else {    
			x[count] = alpha;
			cblas_dscal(k->q[i]-1, alpha/s, &(x[count+1]),1);
		}           
		count += k->q[i];
	}
	/* project onto PSD cone */
	for (i=0; i < k->ssize; ++i){
		projectsdc(&(x[count]),k->s[i],w);
		count += k->s[i]^2;
	}
	/* project onto OTHER cones */
}

void projectsdc(double *X, int n, Work * w)
{ /* project onto the dual positive semi-definite cone */
  int i, j, m=0;
  double * Xs = w->Xs;
  double * Z = w->Z;
  double * e = w->e;
  memcpy(Xs,X,n*n*sizeof(double));

  // Xs = X + X'
  /* save division by 2 until after eigendecomp */
  for (i = 0; i < n; ++i){
    cblas_daxpy(n, 1, &(X[i]), n, &(Xs[i*n]), 1);
  }
  
  // MOVE THIS:
  double EIG_TOL = 1e-4;
  double vlower = -cblas_dnrm2(n*n, Xs, 1);
  LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'V', 'U',\
      n, Xs, n, vlower, 0.0, -1, -1, EIG_TOL, &m, e, Z, n , NULL);

  /* find negative eigen-decomposition */
  memset(Xs, 0, n*n*sizeof(double));
  for (i = 0; i < m; ++i) {
    //cblas_dger(CblasColMajor,n,n, -e[i]/2, &(Z[i*n]), 1, &(Z[i*n]), 1, Xs, n);
    cblas_dsyr(CblasColMajor, CblasLower, n, -e[i]/2, &(Z[i*n]), 1, Xs, n);
  }
  /* fill in upper half */
  for (i = 0; i < n; ++i){   
    for (j = i+1; j < n; ++j){   
      Xs[i + j*n] = Xs[j + i*n];    
    }   
  }   
  cblas_daxpy(n*n, 1, Xs, 1, X, 1);
}
