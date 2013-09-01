#include "private.h"

// forward declare
int LDLInit(cs * A, int P[], double **info);
int LDLFactor(cs * A, int P[], int Pinv[], cs ** L, double **D);
void LDLSolve(double *x, double b[], cs * L, double D[], int P[], double * bp);
int factorize(Data * d,Work * w);

void freePriv(Work * w){
	cs_spfree(w->p->L);coneOS_free(w->p->P);coneOS_free(w->p->D);coneOS_free(w->p->bp);
	coneOS_free(w->p);
}

void solveLinSys(Data *d, Work * w, double * b, const double * s){
  // returns solution to linear system
  // Ax = b with solution stored in b
  LDLSolve(b, b, w->p->L, w->p->D, w->p->P, w->p->bp);
}

int privateInitWork(Data * d, Work * w){ 
	memcpy(w->method, "sparse-direct", 14);
	int n_plus_m = d->n + d->m;
	
	w->p = coneOS_malloc(sizeof(Priv));
	
	
	w->p->P = coneOS_malloc(sizeof(int)*n_plus_m);
	w->p->L = coneOS_malloc(sizeof (cs));
	w->p->bp = coneOS_malloc(n_plus_m * sizeof(double));
	w->p->L->m = n_plus_m;
	w->p->L->n = n_plus_m;
	w->p->L->nz = -1; 
	
	cholmod_allocate_sparse();

	cholmod_sparse w->A->i =;
	cholmod_sparse w->A->i =;
	cholmod_sparse w->A->i =;

	cholmod_dense *x, *b, *r ;
	cholmod_factor *L ;
	cholmod_common c ;
	cholmod_start (&c) ;

	L = cholmod_analyze (A, &c) ;
	cholmod_factorize_p (A, L, &c) ;
	
	return(factorize(d,w));

	NESDIS;
	cholmod_solve2(CHOLMOD_A,L, b, &c,);
	cholmod_sdmult();
}

void LDLSolve(double *x, double b[], cs * L, double D[], int P[], double * bp)
{
  // solves PLDL'P' x = b for x
  int n = L->n;
  if (P == NULL) {
    if (x != b) // if they're different addresses
      memcpy(x,b, n*sizeof(double)); 
    ldl_lsolve(n, x, L->p, L->i, L->x);
    ldl_dsolve(n, x, D); 
    ldl_ltsolve(n, x, L->p, L->i, L->x);
  } else {
    ldl_perm(n, bp, b, P); 
    ldl_lsolve(n, bp, L->p, L->i, L->x);
    ldl_dsolve(n, bp, D); 
    ldl_ltsolve(n, bp, L->p, L->i, L->x);
    ldl_permt(n, x, bp, P); 
	}   
}
