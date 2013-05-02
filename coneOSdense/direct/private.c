#include "private.h"

void privateInitWork(Data * d, Work * w){
	tic();
	int k,j, n=d->n, m=d->m;
	w->p = coneOS_malloc(sizeof(Priv));
	double * A = d->A;
	/* sparse A input: 
	   w->p->A = coneOS_calloc(m*n, sizeof(double));
	   for (j = 0; j < n; ++j) {                     
	   for (k = d->Ap[j]; k < d->Ap[j+1]; ++k) { 
	   w->p->A[d->Ai[k]+j*m] =  d->Ax[k];
	   }
	   }                  
	 */
	// form Gram matrix I + A'A
	w->p->L = coneOS_calloc(n*n, sizeof(double));
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,n, n, m, 1, A, d->m, A, m, 0, w->p->L, n);
	//b_dgemm('T', 'N', n, n, m, 1, A, d->m, A, m, 0, w->p->L, n);

	for (j = 0; j < n; j++) {
		w->p->L[j*n + j] += 1;
	}  
	LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',d->n,w->p->L,d->n);
	// copy L into top half for faster solve steps
	for (k = 0; k < n; ++k){   
		for (j = k+1; j < n; ++j){   
			w->p->L[k + j*n] = w->p->L[j + k*n];		
		}   
	}   
	//coneOS_printf("Factorization time is %4.8f ms\n", tocq());
}

void freePriv(Work * w){
	coneOS_free(w->p->L);
	coneOS_free(w->p);
}

void solveLinSys(Data *d, Work * w, double * b, const double * s){
	// solves Mx = b, for x but stores result in b
	double * A = d->A;
	double * L = w->p->L;
	cblas_dgemv(CblasColMajor, CblasTrans, d->m, d->n, 1, A, d->m, &(b[d->n]), 1, 1, b, 1);
	//b_dgemv('T', d->m, d->n, 1, A, d->m, &(b[d->n]), 1, 1, b, 1);	

	/* lapack dportrs is slower than cblas dtrsv: */
	//LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', d->n, 1, L, d->n, b, d->n);

	/* Solve using forward-substitution, L c = b */
	cblas_dtrsv (CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, d->n, L, d->n, b, 1); 
	//b_dtrsv('L', 'N', 'N', d->n, L, d->n, b, 1); 
	/* Perform back-substitution, U x = c */
	cblas_dtrsv (CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, d->n, L, d->n, b, 1); 
	//b_dtrsv('U', 'N', 'N', d->n, L, d->n, b, 1); 

	cblas_dgemv(CblasColMajor, CblasNoTrans, d->m, d->n, 1, A, d->m,  b, 1, -1, &(b[d->n]), 1);
	//b_dgemv('N', d->m, d->n, 1, A, d->m,  b, 1, -1, &(b[d->n]), 1);
}
