#include "private.h"

int privateInitWork(Data * d, Work * w){
  int len = 43;
  char str[len]; 
  sprintf(str, "dense-indirect, CG: iters %i, tol %.2e", d->CG_MAX_ITS, d->CG_TOL);
  memcpy(w->method, str, len);	
  int j;
  w->p = coneOS_malloc(sizeof(Priv));
	w->p->x = coneOS_malloc(d->n*sizeof(double));
	w->p->p = coneOS_malloc(d->n*sizeof(double));
	w->p->Ap = coneOS_malloc(d->n*sizeof(double));

  double * A = d->Ax;
	/* sparse A format
  w->p->A = coneOS_calloc(d->m*d->n, sizeof(double));
	int k;
	for (j = 0; j < d->n; j++) {                     
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) { 
			w->p->A[d->Ai[k]+j*d->m] =  d->Ax[k];
		}
	}
  */
	// form Gram matrix I+A'A
	w->p->G = coneOS_malloc(d->n*d->n*sizeof(double));
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,d->n, d->n, d->m, 1, A, d->m, A,d->m, 0, w->p->G,d->n);
  //b_dgemm('T', 'N' ,d->n, d->n, d->m, 1, A, d->m, A,d->m, 0, w->p->G,d->n);
  for (j = 0; j < d->n; j++) 
	{ 
		//w->p->G[j*d->n + j] += 1;
		w->p->G[j*d->n + j] += d->RHO_X;
	}
	return(0);
}

void freePriv(Work * w){
	coneOS_free(w->p->p);
	coneOS_free(w->p->x);
	coneOS_free(w->p->Ap);
	coneOS_free(w->p->G);
	coneOS_free(w->p);
}

static void cgCustom(Data *d, Work *w, const double *s, int max_its, double tol);

void solveLinSys(Data *d, Work * w, double * b, const double * s){
	// solves Mx = b, for x but stores result in b
	// s contains warm-start (if available)
	w->p->r = b;
	double * r = w->p->r;
	double * A = d->Ax;
	double * x = w->p->x;
	cblas_dgemv(CblasColMajor, CblasTrans, d->m, d->n, 1, A, d->m, &(b[d->n]),1, 1, r, 1);	
	//b_dgemv('T', d->m, d->n, 1, A, d->m, &(b[d->n]),1, 1, r, 1);	
	cgCustom(d, w, s, d->CG_MAX_ITS, d->CG_TOL);
	cblas_dgemv(CblasColMajor, CblasNoTrans, d->m, d->n, 1, A, d->m, x,1, -1, &(b[d->n]), 1);
	//b_dgemv('N', d->m, d->n, 1, A, d->m, x,1, -1, &(b[d->n]), 1);
	memcpy(b, x, d->n*sizeof(double));
}


static void cgCustom(Data *d, Work *w, const double * s, int max_its, double tol){
	/* solves (I+A'A)x = b */
	/* warm start cg with s */  
	int i = 0, n = d->n;
	double *x = w->p->x;
	double *p = w->p->p; // cg direction
	double *Ap = w->p->Ap; // updated CG direction
	double *r = w->p->r; // cg residual

	double *G = w->p->G; // Gram matrix

	double alpha, beta, rsnew=0;
	if (s==NULL){
		memset(x,0,n*sizeof(double));
	}
	else{
		memcpy(x,s,n*sizeof(double));
		cblas_dsymv(CblasColMajor, CblasUpper,n, -1, G, n, x,1, 1, r, 1);
		//b_dsymv('U', n, -1, G, n, x,1, 1, r, 1);
	}
	memcpy(p, r, n*sizeof(double));
	//double rsold=cblas_dnrm2(n,r,1);
	double rsold=calcNorm(r,n);
	for (i=0; i< max_its; i++){
		cblas_dsymv(CblasColMajor, CblasUpper,n, 1, G, n, p, 1, 0, Ap, 1);
		//b_dsymv('U', n, 1, G, n, p, 1, 0, Ap, 1);

		//beta = cblas_ddot(n, p, 1, Ap, 1);
		beta = innerProd(p,Ap,n);
		alpha=(rsold*rsold)/beta;

		addScaledArray(x,p,n,alpha);
		//cblas_daxpy(n,alpha,p,1,x,1);
		addScaledArray(r,Ap,n,-alpha);
		//cblas_daxpy(n,-alpha,Ap,1,r,1);

		//rsnew=cblas_dnrm2(n,r,1);
		rsnew=calcNorm(r,n);
		if (rsnew<tol){
			break;
		}
		scaleArray(p,(rsnew*rsnew)/(rsold*rsold),n);
		//cblas_dscal(n,(rsnew*rsnew)/(rsold*rsold),p,1);
		addScaledArray(p,r,n,1);
		//cblas_daxpy(n,1,r,1,p,1);
		rsold=rsnew;
	} 
	//printf("terminating cg residual = %4f, took %i itns\n",rsnew,i);
}
