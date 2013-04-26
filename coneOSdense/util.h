#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "coneOS.h"
#include "cones.h"

void tic(void); 
float toc(void); 
float tocq(void); 
void printConeData(Cone * k);
void printData(Data * d);
void printAll(Data * d, Work * w);


#endif
