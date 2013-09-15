#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#include "coneOS.h"
#include <cblas.h>
#include <math.h>

void setAsScaledArray(double *x, const double * a,const double b,int len);
void scaleArray(double * a,const double b,int len);
double innerProd(const double * x, const double * y, int len);
double calcNorm(const double * v,int len);
double calcNormInf(const double *a, int l);
double calcNormSq(const double * v,int len);
void addScaledArray(double * a, const double * b, int n, const double sc);
double calcNormDiff(const double *a,const double *b, int l);
double calcNormInfDiff(const double *a, const double *b, int l);
void accumByAtrans(Data * d, const double *x, double *y);
void accumByA(Data * d, const double *x, double *y);
#endif
