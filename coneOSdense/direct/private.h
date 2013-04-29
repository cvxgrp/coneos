#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "coneOS.h"
#include "linAlg.h"
#include <lapacke.h>

struct PRIVATE_DATA{
	double * L;
};

//Work * initWork(Data* d);
//void formQ(Data * d, Work * w);
//void freePriv(Work * w);
//void projectLinSys(Data * d,Work * w);
//void cgCustom(cs *Q,double * b,double * x,int max_its,double tol);
//void multByQ(const cs *Q, const double *x, double *y);

#endif
