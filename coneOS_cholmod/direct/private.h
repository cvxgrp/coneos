#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "cholmod.h"
#include "coneOS.h"

struct PRIVATE_DATA {
	cs * L; /* KKT, and factorization matrix L resp. */
	double * D; /* diagonal matrix of factorization */
	int * P; /* permutation of KKT matrix for factorization */
	double * bp; /* workspace memory for solves */
};

extern void amd_info (double Info [ ]);
extern int amd_order (int n, const int Ap [ ], const int Ai [ ], int P [ ], double Control [ ], double Info [ ]);

#endif
