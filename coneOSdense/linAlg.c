#include "linAlg.h"

// x = b*a
inline void setAsScaledArray(double *x, const double * a,const double b,int len);

// a*= b
inline void scaleArray(double * a,const double b,int len);

// x'*y
inline double innerProd(const double * x, const double * y, int len);

inline double calcNormDiff(double *a, double *b, int l); 

// ||v||_2^2
inline double calcNormSq(const double * v,int len);

// ||v||_2
inline double calcNorm(const double * v,int len);

// daxpy
inline void addScaledArray(double * a, const double * b, int n, const double sc);

static inline void accumByAtrans(Data * d, const double *x, double *y); 

static inline void accumByA(Data * d, const double *x, double *y); 
