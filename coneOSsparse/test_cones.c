#include <stdio.h>
#include <stdlib.h>
#include "cones.h"
#include "coneOS.h"
#include "linAlg.h"

int main(int argc, char **argv)
{
    double v[3] = {-100.4138, 129.3833, -66.3384};
    printf("clear all; i=1;\n");
    for(int i=0; i<1e4; i++) {
        v[0] = (double)rand()/(double)RAND_MAX - (double)rand()/(double)RAND_MAX;
        v[1] = (double)rand()/(double)RAND_MAX - (double)rand()/(double)RAND_MAX;
        v[2] = (double)rand()/(double)RAND_MAX - (double)rand()/(double)RAND_MAX;
        printf("p(i,:) = [%4f, %4f, %4f];\n",v[0],v[1],v[2]);
        projExpCone(v);
        printf("po(i,:) = [%4f, %4f, %4f];\n",v[0],v[1],v[2]);
        printf("px(i,:) = project_exp_bisection(p(i,:));\n");
        printf("nm(i) = norm(px(i,:) -po(i,:));\n");
        printf("i = i+1;\n");
    }
}
