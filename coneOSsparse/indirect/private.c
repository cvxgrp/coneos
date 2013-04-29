#include "private.h"
#include "omp.h"

static inline void calcAx(Data * d, Work * w, const double * x, double * y);
static void cgCustom(Data *d, Work *w, const double *s, double * b, int max_its, double tol);
static inline void accumByA(Data * d, Work * w, const double *x, double *y);
static inline void accumByAtrans(Data *d, Work * w, const double *x, double *y);
static inline void _accumAtrans(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y);
static inline void transpose (Data * d, Work * w);


void privateInitWork(Data * d, Work * w){
  w->p = coneOS_malloc(sizeof(Priv));
  w->p->p = coneOS_malloc((d->n)*sizeof(double));
  w->p->r = coneOS_malloc((d->n)*sizeof(double));
  w->p->Ap = coneOS_malloc((d->n)*sizeof(double));
  w->p->tmp = coneOS_malloc((d->m)*sizeof(double));

  w->p->Ati = coneOS_malloc((d->Ap[d->n])*sizeof(int));
  w->p->Atp = coneOS_malloc((d->m+1)*sizeof(int));
  w->p->Atx = coneOS_malloc((d->Ap[d->n])*sizeof(double));
  transpose(d,w);
}

static inline void transpose (Data * d, Work * w)
{
  int * Ci = w->p->Ati;
  int * Cp = w->p->Atp;
  double * Cx = w->p->Atx;
  int m = d->m;
  int n = d->n;

  int * Ap = d->Ap;
  int * Ai = d->Ai;
  double * Ax = d->Ax;

  int p, j, q, *z;
  z = coneOS_calloc(m,sizeof(int));
  for (p = 0 ; p < Ap [n] ; p++) z [Ai [p]]++ ;          /* row counts */
  cs_cumsum (Cp, z, m) ;                                 /* row pointers */
  for (j = 0 ; j < n ; j++)
  {
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      Ci [q = z [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
      if (Cx) Cx [q] = Ax [p] ;
    }
  }
  coneOS_free(z);
}

void freePriv(Work * w){
	coneOS_free(w->p->p);
	coneOS_free(w->p->r);
	coneOS_free(w->p->Ap);
	coneOS_free(w->p->tmp);
  coneOS_free(w->p->Ati);
  coneOS_free(w->p->Atx);
  coneOS_free(w->p->Atp);
	coneOS_free(w->p);
}

void solveLinSys(Data *d, Work * w, double * b, const double * s){
	// solves Mx = b, for x but stores result in b
	// s contains warm-start (if available)
	accumByAtrans(d,w, &(b[d->n]), b);
	// solves (I+A'A)x = b, s warm start, solution stored in b
	cgCustom(d, w, s, b, d->CG_MAX_ITS, d->CG_TOL);
	scaleArray(&(b[d->n]),-1,d->m);
	accumByA(d, w, b, &(b[d->n]));
}

static void cgCustom(Data *d, Work *w, const double * s, double * b, int max_its, double tol){
	/* solves (I+A'A)x = b */
	/* warm start cg with x */  
	int i = 0, n = d->n;
	double *p = w->p->p; // cg direction
	double *Ap = w->p->Ap; // updated CG direction
	double *r = w->p->r; // cg residual

	double alpha, rsnew=0;
	if (s==NULL){
		memcpy(r,b,n*sizeof(double));
		memset(b,0.0,n*sizeof(double));
	}
	else{
		calcAx(d,w,s,r);
		addScaledArray(r,b,n,-1);
		scaleArray(r,-1,n);
		memcpy(b,s,n*sizeof(double));
	}
	memcpy(p,r,n*sizeof(double));
	double rsold=calcNorm(r,n);

	for (i=0; i< max_its; i++){
		calcAx(d,w,p,Ap);
		
		alpha=(rsold*rsold)/innerProd(p,Ap,n);
		addScaledArray(b,p,n,alpha);
		addScaledArray(r,Ap,n,-alpha);   

		rsnew=calcNorm(r,n);
		if (rsnew<tol){
			break;
		}   
		scaleArray(p,(rsnew*rsnew)/(rsold*rsold),n);
		addScaledArray(p,r,n,1);
		rsold=rsnew;
	} 
	//printf("terminating cg residual = %4f, took %i itns\n",rsnew,i);
}

static inline void calcAx(Data * d, Work * w, const double * x, double * y){
	double * tmp = w->p->tmp;
	memset(tmp,0,d->m*sizeof(double));
	accumByA(d,w,x,tmp);
	memset(y,0,d->n*sizeof(double));
	accumByAtrans(d,w,tmp,y);
	addScaledArray(y,x,d->n,1);
}

static inline void accumByA(Data * d, Work * w, const double *x, double *y)
{
  _accumAtrans(d->m,w->p->Atx,w->p->Ati,w->p->Atp,x,y);
}
static inline void accumByAtrans(Data *d, Work * w, const double *x, double *y){
  _accumAtrans(d->n,d->Ax,d->Ai,d->Ap,x,y);
}

static inline void _accumAtrans(int n, double * Ax, int * Ai, int * Ap, const double *x, double *y)
{
	/* y  = A'*x 
	A in column compressed format 
	parallelizes over columns (rows of A')
	*/
	int p, j;
	int c1, c2;
	double yj;
#pragma omp parallel for private(p,c1,c2,yj) 
	for (j = 0 ; j < n ; j++)
	{
		yj = y[j];
		c1 = Ap[j]; c2 = Ap[j+1];
		for (p = c1 ; p < c2 ; p++)        
		{   
			yj += Ax[p] * x[ Ai[p] ] ;
		}
		y[j] = yj;
	}
}

//
//void accumByA(Data * d, const double *x, double *y)
//{
//*
///*y  = A*x 
//	A in column compressed format  
//	this parallelizes over columns and uses
//	pragma atomic to prevent concurrent writes to y 
//	 */
//	int p, j, n, *Ap, *Ai ;
//	double *Ax ;
//	n = d->n ; Ap = d->Ap ; Ai = d->Ai ; Ax = d->Ax ;
//	int c1, c2;
//	double xj;
////#pragma omp parallel for private(p,c1,c2,xj) 
//	for (j = 0 ; j < n ; j++)
//	{
//		xj = x[j];
//		c1 = Ap[j]; c2 = Ap[j+1];
//		for (p = c1 ; p < c2 ; p++)        
//		{
////#pragma omp atomic
//			y [Ai[p]] += Ax [p] * xj ;
//		}
//	}
//}
//
//void accumByAtrans(Data * d, const double *x, double *y)
//{
//	/* y  = A'*x 
//	A in column compressed format 
//	parallelizes over columns (rows of A')
//	*/
//	int p, j, n, *Ap, *Ai ;
//	double *Ax ;
//	n = d->n ; Ap = d->Ap ; Ai = d->Ai ; Ax = d->Ax ;
//	/* parallel matrix multiply */
//	int c1, c2;
//	double yj;
//#pragma omp parallel for private(p,c1,c2,yj) 
//	for (j = 0 ; j < n ; j++)
//	{
//		yj = y[j];
//		c1 = Ap[j]; c2 = Ap[j+1];
//		for (p = c1 ; p < c2 ; p++)        
//		{   
//			yj += Ax[p] * x[ Ai[p] ] ;
//		}
//		y[j] = yj;
//	}
//}
