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
  //b_dscal(len,b,a,1);
  cblas_dscal(len,b,a,1);
}

// x'*y
inline double innerProd(const double * x, const double * y, int len){
  //return b_ddot(len, x,1,y,1);
  return cblas_ddot(len, x,1,y,1);
}

// ||v||_2
inline double calcNorm(const double * v,int len){
  //return b_dnrm2(len, v, 1);
  return cblas_dnrm2(len, v, 1);
}

inline double calcNormInf(const double *a, int l){
	double tmp, max = 0.0;
	int i;
	for ( i=0; i<l; ++i){
		tmp = fabs(a[i]);
		if(tmp > max) max = tmp;
	}
	return max;
}

// ||v||_2^2
inline double calcNormSq(const double * v,int len){
  double nrm = calcNorm(v,len);
  return nrm*nrm;
}

// daxpy a += sc*b
inline void addScaledArray(double * a, const double * b, int n, const double sc){
  //b_daxpy(n, sc, b, 1, a, 1);
  cblas_daxpy(n, sc, b, 1, a, 1);
}

inline double calcNormDiff(const double *a,const double *b, int l) {
  double nmDiff = 0.0, tmp;
  int i;
  for ( i=0; i<l; ++i ){
    tmp = (a[i] - b[i]);
    nmDiff += tmp * tmp;
  }  
  return sqrt(nmDiff);
}

inline double calcNormInfDiff(const double *a, const double *b, int l) {
	double tmp, max = 0.0;
	int i;
	for ( i=0; i<l; ++i){
		tmp = fabs(a[i] - b[i]);
		if(tmp > max) max = tmp;
	}
	return max;
}

static inline void accumByAtrans(Data * d, const double *x, double *y) 
{
	cblas_dgemv(CblasColMajor, CblasTrans, d->m, d->n, 1.0, d->Ax, d->m, x, 1, 1.0, y, 1);
}

static inline void accumByA(Data * d, const double *x, double *y) 
{
 	cblas_dgemv(CblasColMajor, CblasNoTrans, d->m, d->n, 1.0, d->Ax, d->m, x, 1, 1.0, y, 1);
}

/*

inline void b_daxpy(const int n, const double alpha, const double *x, const int incx, double *y, const int incy)
{
  daxpy(&n,&alpha,x,&incx,y,&incy);
  //cblas_daxpy(n,alpha,x,incx,y,incy);
}

inline void b_dsyr(char Uplo, const int N, const double alpha, const double *X, const int incX, double *A, const int lda)
{
  //dsyr(&Uplo,&N,&alpha,X,&incX,A,&lda);
  cblas_dsyr(CblasColMajor, Uplo, N, alpha, X, incX, A, lda);
}

inline void b_dscal(const int N, const double alpha, double *X, const int incX)
{
  //dscal(&N,&alpha,X,&incX);
  cblas_dscal(N,alpha,X,incX);
}

inline double b_ddot(const int n, const double *x, const int incx, const double *y, const int incy)
{
  //return ddot(&n,x,&incx,y,&incy);
  return cblas_ddot(n,x,incx,y,incy);
}

inline double b_dnrm2(const int N, const double *X, const int incX)
{
  return cblas_dnrm2(N,X,incX);
  //return dnrm2(&N,X,&incX);
}

inline void b_dgemm(
              char TransA,
              char TransB, 
              const int M, 
              const int N,
              const int K,
              const double alpha, 
              const double *A, 
              const int lda, 
              const double *B, 
              const int ldb, 
              const double beta, 
              double *C, 
              const int ldc){
  dgemm(&TransA, &TransB, &M, &N, &K,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
  //cblas_dgemm();
}

inline void b_dgemv(char Trans, 
                    const int m,
                    const int n,
                    const double alpha, 
                    const double  *a, 
                    const int lda,  
                    const double  *x,
                    const int incx,  
                    const double beta,  
                    double  *y, 
                    const int incy){
  dgemv(&Trans,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}

inline void b_dtrsv(char Uplo, char TransA, char Diag, const int N, const double *A, const int lda, double *X, const int incX)
{
  dtrsv(&Uplo, &TransA, &Diag, &N, A, &lda, X, &incX);
}

inline void b_dsymv(char Uplo, const int N, const double alpha, const double *A, 
                    const int lda, const double *X, const int incX, const double beta,
                    double *Y, const int incY){
  dsymv(&Uplo, &N, &alpha, A, &lda, X, &incX, beta, Y, &incY);
}

*/

#endif
