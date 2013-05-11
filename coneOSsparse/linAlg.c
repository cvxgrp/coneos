#include "linAlg.h"
/*
 * All basic linear operations are inlined (and further optimized) by the 
 * compiler. If compiling without optimization, causes code bloat.
 */

// x = b*a
inline void setAsScaledArray(double *x, const double * a,const double b,int len); 

// a*= b
inline void scaleArray(double * a,const double b,int len);

// x'*y
inline double innerProd(const double * x, const double * y, int len);

// ||v||_2^2
inline double calcNormSq(const double * v,int len);

// ||v||_2
inline double calcNorm(const double * v,int len);

// saxpy a += sc*b
inline void addScaledArray(double * a, const double * b, int n, const double sc);

inline double calcNormDiff(double *a, double *b, int l);

static inline void accumByAtrans(Data * d, const double *x, double *y); 

static inline void accumByA(Data * d, const double *x, double *y);
