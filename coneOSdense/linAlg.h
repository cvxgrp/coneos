#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#include "coneOS.h"
#include <cblas.h>
#include <math.h>

// x = b*a
inline void setAsScaledArray(double *x, const double * a,const double b,int len) {
  int i;
  for( i=0;i<len;++i ) x[i] = b*a[i];
}

// a*= b
inline void scaleArray(double * a,const double b,int len){
  cblas_dscal(len,b,a,1);
}

// x'*y
inline double innerProd(const double * x, const double * y, int len){
  return cblas_ddot(len, x,1,y,1);
}

// ||v||_2
inline double calcNorm(const double * v,int len){
  return cblas_dnrm2(len, v, 1);
}


// ||v||_2^2
inline double calcNormSq(const double * v,int len){
  double nrm = calcNorm(v,len);
  return nrm*nrm;
}

// daxpy a += sc*b
inline void addScaledArray(double * a, const double * b, int n, const double sc){
  cblas_daxpy(n, sc, b, 1, a, 1);
}

inline double calcNormDiff(double *a, double *b, int l) {
    double nmDiff = 0.0, tmp;
    int i;
    for ( i=0; i<l; ++i ){
        tmp = (a[i] - b[i]);
		nmDiff += tmp * tmp;
	}  
  return sqrt(nmDiff);
}

#endif
