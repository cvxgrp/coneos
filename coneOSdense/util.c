#include "util.h"

static clock_t tic_timestart;

void tic(void) {
	tic_timestart = clock();
}

float toc(void) {
	clock_t tic_timestop;
	tic_timestop = clock();
	coneOS_printf("time: %8.4f seconds.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
	return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

float tocq(void) {
	clock_t tic_timestop;
	tic_timestop = clock();
	return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

void printConeData(Cone * k){
	int i;
	coneOS_printf("num zeros = %i\n",k->f);
	coneOS_printf("num LP = %i\n",k->l);
	coneOS_printf("num SOCs = %i\n",k->qsize);
	coneOS_printf("soc array:\n");
	for ( i=0;i<k->qsize;i++){
		coneOS_printf("%i\n",k->q[i]);
	}
	/*
	coneOS_printf("num SDCs = %i\n",k->ssize);
	coneOS_printf("sdc array:\n");
	for ( i=0;i<k->ssize;i++){
		coneOS_printf("%i\n",k->s[i]);
	}
	*/
}

void printData(Data * d){
	coneOS_printf("d->n is %i\n",d->n);
	coneOS_printf("d->m is %i\n",d->m);
	coneOS_printf("d->b[0] is %4f\n",d->b[0]);
	coneOS_printf("d->c[0] is %4f\n",d->c[0]);
	coneOS_printf("d->MAX_ITERS is %i\n",d->MAX_ITERS);
	coneOS_printf("d->CG_MAX_ITS is %i\n",d->CG_MAX_ITS);
	coneOS_printf("d->ALPH is %6f\n",d->ALPH);
	coneOS_printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
	//coneOS_printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
	//coneOS_printf("d->UNDET_TOL is %6f\n",d->UNDET_TOL);
	coneOS_printf("d->CG_TOL is %6f\n",d->CG_TOL);
	coneOS_printf("d->A[0,0] is %4f\n",d->A[0]);
}

void printAll(Data * d, Work * w){
	int i;
	coneOS_printf("\n u_t is \n");
	for( i=0;i<w->l;i++){
		coneOS_printf("%f\n",w->u_t[i]);
	}
	coneOS_printf("\n u is \n");
	for( i=0;i<w->l;i++){
		coneOS_printf("%f\n",w->u[i]);
	}
	coneOS_printf("\n v is \n");
	for( i=0;i<w->l;i++){
		coneOS_printf("%f\n",w->v[i]);
	}
}

