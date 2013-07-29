#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "cs.h"
#include "amd.h"
#include "ldl.h"
#include "coneOS.h"

struct PRIVATE_DATA {
	cs * L; /* KKT, and factorization matrix L resp. */
	double * D; /* diagonal matrix of factorization */
	int * P; /* permutation of KKT matrix for factorization */
};

// XXX: should be named LDL
// also, these routines don't need to be "public"
//void LDLInit(cs * A, int P[], double **info);
//void LDLFactor(cs * A, int P[], int Pinv[], cs ** L, double **D);
//void LDLSolve(double *x, double b[], cs * L, double D[], int P[]);
//int privateInitWork(Data * d, Work * w);
//void solveLinSys(Data * d, Work * w, double * b, double * w);
//void freePriv(Work * w);
//cs * formKKT(Data * d, Work * w);
//void factorize(Data * d,Work * w);

extern void amd_info (double Info [ ]);
extern int amd_order (int n, const int Ap [ ], const int Ai [ ], int P [ ], double Control [ ], double Info [ ]);

#endif
