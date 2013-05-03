#ifndef CONES_H_GUARD                                                              
#define CONES_H_GUARD

#include <string.h>

typedef struct Cone_t {
    int f;          /* number of linear equality constraints */
    int l;          /* length of LP cone */
    int *q;   		/* array of second-order cone constraints */
    int qsize;      /* length of SOC array */
	int *s;			/* array of SD constraints */
	int ssize;		/* length of SD array */
} Cone;

#include "coneOS.h"
void projCone(double *x,Cone *k, Work * w);
#endif
