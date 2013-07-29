#include "private.h"

// forward declare
int LDLInit(cs * A, int P[], double **info);
int LDLFactor(cs * A, int P[], int Pinv[], cs ** L, double **D);
void LDLSolve(double *x, double b[], cs * L, double D[], int P[]);
int factorize(Data * d,Work * w);

void freePriv(Work * w){
	cs_spfree(w->p->L);coneOS_free(w->p->P);coneOS_free(w->p->D);
	coneOS_free(w->p);
}

void solveLinSys(Data *d, Work * w, double * b, const double * s){
  // returns solution to linear system
  // Ax = b with solution stored in b
  LDLSolve(b, b, w->p->L, w->p->D, w->p->P);
}

int privateInitWork(Data * d, Work * w){ 
	memcpy(w->method, "sparse-direct", 14);
	int n_plus_m = d->n + d->m;
	w->p = coneOS_malloc(sizeof(Priv));
	w->p->P = coneOS_malloc(sizeof(int)*n_plus_m);
	w->p->L = coneOS_malloc(sizeof (cs));
	w->p->L->m = n_plus_m;
	w->p->L->n = n_plus_m;
	w->p->L->nz = -1; 
	return(factorize(d,w));
}

cs * formKKT(Data * d, Work * w){
	/* ONLY UPPER TRIANGULAR PART IS STUFFED
	 * forms column compressed KKT matrix
	 * assumes column compressed form A matrix
	 *
	 * forms upper triangular part of [I A'; A -I]
	 */
	int j, k, kk;
	/* I at top left */
	const int Anz = d->Ap[d->n];
	const int Knzmax = d->n + d->m + Anz;
	cs * K = cs_spalloc(d->m + d->n, d->m + d->n, Knzmax, 1, 1);
	kk = 0;
	for (k = 0; k < d->n; k++){
		K->i[kk] = k;
		K->p[kk] = k;
		K->x[kk] = d->RHO_X;
		kk++;
	}
	/* A^T at top right : CCS: */
	for (j = 0; j < d->n; j++) {                 
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) { 
			K->p[kk] = d->Ai[k] + d->n;
			K->i[kk] = j;
			K->x[kk] = d->Ax[k];
			kk++;
		}   
	}
	/* -I at bottom right */
	for (k = 0; k < d->m; k++){
		K->i[kk] = k + d->n;
		K->p[kk] = k + d->n;
		K->x[kk] = -1;
		kk++;
	}
	// assert kk == Knzmax
	K->nz = Knzmax;
	cs * K_cs = cs_compress(K);
	cs_spfree(K);
	return(K_cs);
}


int factorize(Data * d,Work * w){
	//tic();
	cs * K = formKKT(d,w);
	//if(d->VERBOSE) coneOS_printf("KKT matrix factorization info:\n");
	double *info;
	int amd_status = LDLInit(K, w->p->P, &info);
	if (amd_status < 0) return(amd_status);
	/*
  if(d->VERBOSE) {
#ifdef DLONG
		amd_l_info(info);
#else
		amd_info(info);
#endif
	}
  */
	int * Pinv = cs_pinv(w->p->P, w->l-1);
	cs * C = cs_symperm(K, Pinv, 1); 
	int ldl_status = LDLFactor(C, NULL, NULL, &w->p->L, &w->p->D);
	//if(d->VERBOSE) coneOS_printf("KKT matrix factorization took %4.8f ms\n",tocq());
	cs_spfree(C);cs_spfree(K);coneOS_free(Pinv);coneOS_free(info);
	return(ldl_status);
}

int LDLInit(cs * A, int P[], double **info) {
	*info  = (double *) coneOS_malloc(AMD_INFO * sizeof(double));
#ifdef DLONG
	return(amd_l_order(A->n, A->p, A->i, P, (double *) NULL, *info));
#else
	return(amd_order(A->n, A->p, A->i, P, (double *) NULL, *info));
#endif
}

int LDLFactor(cs * A, int P[], int Pinv[], cs **L , double **D) 
{
	int n = A->n;
	(*L)->p = (int *) coneOS_malloc((1 + n) * sizeof(int));
	
	//int Parent[n], Lnz[n], Flag[n], Pattern[n];
	//double Y[n];

	int * Parent = coneOS_malloc(n*sizeof(int));
	int * Lnz = coneOS_malloc(n*sizeof(int));
	int * Flag = coneOS_malloc(n*sizeof(int));
	int * Pattern = coneOS_malloc(n*sizeof(int));
	double * Y = coneOS_malloc(n*sizeof(double));

	ldl_symbolic(n, A->p, A->i, (*L)->p, Parent, Lnz, Flag, P, Pinv);

	(*L)->nzmax = *((*L)->p + n);
	(*L)->x = (double *) coneOS_malloc((*L)->nzmax * sizeof(double));
	(*L)->i =    (int *) coneOS_malloc((*L)->nzmax * sizeof(int));
	*D  = (double *) coneOS_malloc(n * sizeof(double));
	
	if(!(*D) || !(*L)->i || !(*L)->x || !Y || !Pattern || !Flag || !Lnz || !Parent) return -1;
	
	int kk = ldl_numeric(n, A->p, A->i, A->x, (*L)->p, Parent, Lnz, (*L)->i, (*L)->x, *D, Y, Pattern, Flag, P, Pinv);

	coneOS_free(Parent);
	coneOS_free(Lnz);
	coneOS_free(Flag);
	coneOS_free(Pattern);
	coneOS_free(Y);
	return(n - kk);
}

void LDLSolve(double *x, double b[], cs * L, double D[], int P[])
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
    double * bp = coneOS_malloc(n * sizeof(double));
    ldl_perm(n, bp, b, P); 
    ldl_lsolve(n, bp, L->p, L->i, L->x);
    ldl_dsolve(n, bp, D); 
    ldl_ltsolve(n, bp, L->p, L->i, L->x);
    ldl_permt(n, x, bp, P); 
	coneOS_free(bp);
  }   
}
