#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#include <math.h>
#include "coneOS.h"

/*
 * All basic linear operations are inlined (and further optimized) by the 
 * compiler. If compiling without optimization, causes code bloat.
 */

// x = b*a
inline void setAsScaledArray(double *x, const double * a,const double b,int len) {
  int i;
  for( i=0;i<len;++i ) x[i] = b*a[i];
}

// a*= b
inline void scaleArray(double * a,const double b,int len){
  int i;
  for( i=0;i<len;++i) a[i]*=b;
}

// x'*y
inline double innerProd(const double * x, const double * y, int len){
  int i;
  double ip = 0.0;
  for ( i=0;i<len;++i){
    ip += x[i]*y[i];
  }
  return ip;
}

// ||v||_2^2
inline double calcNormSq(const double * v,int len){
  int i;
  double nmsq = 0.0;
  for ( i=0;i<len;++i){
    nmsq += v[i]*v[i];
  }
  return nmsq;
}

// ||v||_2
inline double calcNorm(const double * v,int len){
  return sqrt(calcNormSq(v, len));
}

// saxpy a += sc*b
inline void addScaledArray(double * a, const double * b, int n, const double sc){
  int i;
  for (i=0;i<n;++i){
    a[i] += sc*b[i];
  }
}

inline double calcNormDiff(double *a, double *b, int l) {
    double nmDiff = 0.0, tmp;
    int i;
    for ( i=0; i<l; ++i){
        tmp = (a[i] - b[i]);
		nmDiff += tmp * tmp;
	}  
    return sqrt(nmDiff);
}

#endif
