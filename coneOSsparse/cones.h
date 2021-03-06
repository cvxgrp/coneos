#ifndef CONES_H_GUARD                                                              
#define CONES_H_GUARD

#include <string.h>
#include <math.h>

typedef struct Cone_t {
    int f;          /* number of linear equality constraints */
    int l;          /* length of LP cone */
    int *q;   		/* array of second-order cone constraints */
    int qsize;      /* length of SOC array */
	int *s;			/* array of SD constraints */
	int ssize;		/* length of SD array */
    int ep;         /* number of triples in exponential cone */
    int ed;         /* number of triples in dual exponential cone */
} Cone;

#include "coneOS.h"
void projCone(double *x,Cone *k, Work * w, int iter);
void projExpCone(double * v);
#endif
