#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "coneOS.h"
#include <cblas.h>

struct PRIVATE_DATA{
	double * x;
	double * p;
	double * r;  
	double * Ap;
  /* Gram matrix */
	double * G;
};

//Work * initWork(Data* d);
//void formQ(Data * d, Work * w);
//void freePriv(Work * w);
//void projectLinSys(Data * d,Work * w);
//void cgCustom(cs *Q,double * b,double * x,int max_its,double tol);
//void multByQ(const cs *Q, const double *x, double *y);

#endif
