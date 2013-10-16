#ifndef coneOS_H_GUARD                                                              
#define coneOS_H_GUARD

// redefine printfs and memory allocators as needed
#ifdef MATLAB_MEX_FILE
  #include "mex.h"
  #define coneOS_printf   mexPrintf
  #define coneOS_free     mxFree
  #define coneOS_malloc   mxMalloc
  #define coneOS_calloc   mxCalloc
#elif defined PYTHON
  #include <Python.h>
  #include <stdlib.h>
  #define coneOS_printf   PySys_WriteStdout
  #define coneOS_free     free
  #define coneOS_malloc   malloc
  #define coneOS_calloc   calloc
#else
  #include <stdio.h>
  #include <stdlib.h>
  #define coneOS_printf   printf
  #define coneOS_free     free
  #define coneOS_malloc   malloc
  #define coneOS_calloc   calloc
#endif


/* struct that containing standard problem data */
typedef struct PROBLEM_DATA {
  int n, m; /* problem dimensions */
  /* problem data, A, b, c: */
  double * Ax;
  int * Ai, * Ap;
  int Anz;
  double * b, * c;
  int MAX_ITERS, CG_MAX_ITS;
  double EPS_ABS, ALPH, CG_TOL, UNDET_TOL, RHO_X;
  int VERBOSE, NORMALIZE;  // boolean
} Data;

/* contains primal-dual solution vectors */
typedef struct SOL_VARS {
  double * x, * y, *s;
} Sol;

/* contains terminating information */
typedef struct INFO {
	int iter;
	char status[16];
	double pobj;
	double dobj;
	double presid;
	double dresid;
	double gap;
	double time;
} Info;

typedef struct PRIVATE_DATA Priv;

typedef struct WORK {
  double *u, *v, *u_t, *u_prev;
  double *h, *g;  
  double gTh, sc_b, sc_c, scale;
  double nm_b, nm_c, nm_Q;
  double *D, *E;
  Priv * p;
  /* workspace for eigenvector decompositions: */
  double * Xs, *Z, *e;
  int l;
  char method[16];
} Work;

// to hold residual information
struct residuals {
	double resDual;
	double resPri;
    double relGap;
    double tau;
	double kap;
};

#include <string.h>    
#include <sys/time.h>
#include <math.h>
#include "cones.h"
#include "linAlg.h"
#include "util.h"

// these are actually library "api"'s
int coneOS(Data * d, Cone * k, Sol * sol, Info * info);

// these are pulled in from private.o
int privateInitWork(Data * d, Work * w);
// solves [I A';A -I] x = b, stores result in b, s contains warm-start
void solveLinSys(Data * d, Work * w, double * b, const double * s);
void freePriv(Work * w);

#endif
