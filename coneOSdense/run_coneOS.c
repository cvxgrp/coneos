#include "coneOS.h"
#include "run_coneOS.h"

#ifndef DEMO_PATH
#define DEMO_PATH "../data_coneOS"
#endif 

#define NUM_TRIALS 1 

int main(int argc, char **argv)
{
	FILE * fp;
	if(open_file(argc, argv, 1, DEMO_PATH, &fp)==-1) return -1;
	Cone * k = malloc(sizeof(Cone));
	Data * d = malloc(sizeof(Data));
	read_in_data(fp,d,k);
	fclose(fp);
	int i;
	Sol * sol = malloc(sizeof(Sol));
	Info * info = malloc(sizeof(Info));
	for (i=0;i<NUM_TRIALS;i++)
	{
		coneOS(d,k,sol,info);
	}
	freeData(d,k);
	freeSol(sol);
	coneOS_free(info);
	return 0;
}

void read_in_data(FILE * fp,Data * d, Cone * k){
	/* MATRIX IN DATA FILE MUST BE DENSE COLUMN MAJOR FORMAT */
	fscanf(fp, "%i", &(d->n));
	fscanf(fp, "%i", &(d->m));
	fscanf(fp, "%i", &(k->f));
	fscanf(fp, "%i", &(k->l));
	fscanf(fp, "%i", &(k->qsize)); 
	fscanf(fp, "%i", &(k->ssize)); 

	k->q = malloc(sizeof(int)*k->qsize);
	for(int i = 0; i < k->qsize; i++)
	{ 
		fscanf(fp, "%i", &k->q[i]);
	}
	k->s = malloc(sizeof(int)*k->ssize);
	for(int i = 0; i < k->ssize; i++)
	{ 
		fscanf(fp, "%i", &k->s[i]);
	}
	d->b = malloc(sizeof(double)*d->m);
	for(int i = 0; i < d->m; i++)
	{ 
		fscanf(fp, "%lf", &d->b[i]);
	}
	d->c = malloc(sizeof(double)*d->n);
	for(int i = 0; i < d->n; i++)
	{ 
		fscanf(fp, "%lf", &d->c[i]);
	}
	fscanf(fp, "%i", &(d->MAX_ITERS));
	fscanf(fp, "%i", &(d->CG_MAX_ITS)); 

	fscanf(fp, "%lf", &(d->ALPH));
	fscanf(fp, "%lf", &(d->UNDET_TOL)); 
	fscanf(fp, "%lf", &(d->EPS_ABS)); 
	fscanf(fp, "%lf", &(d->EPS_REL)); 
	fscanf(fp, "%lf", &(d->CG_TOL));
	fscanf(fp, "%i", &(d->VERBOSE));
	fscanf(fp, "%i", &(d->NORMALIZE));

	d->Anz = d->n*d->m;
	d->Ax = malloc(sizeof(double)*d->Anz);
	for(int i = 0; i < d->Anz; i++)
	{
		fscanf(fp, "%lf", &d->Ax[i]);
	}
	d->Ai = NULL;
	d->Ap = NULL;

	//		fscanf(fp, "%zu", &NNZ);
	//		int *Kr = malloc(sizeof(int)*NNZ);
	//		for(int i = 0; i < NNZ; i++)
	//		{
	//		fscanf(fp, "%i", &Kr[i]);
	//		}
	//		int *Kp=malloc(sizeof(int)*(w->l+1));
	//		for(int i = 0; i < w->l+1; i++)
	//		{
	//		fscanf(fp, "%i", &Kp[i]);
	//		}
	//		double *Kx=malloc(sizeof(double)*NNZ);
	//		for(int i = 0; i < NNZ; i++)
	//		{
	//		fscanf(fp, "%lf", &Kx[i]);
	//		}

}

int open_file(int argc, char ** argv, int idx, char * default_file, FILE ** fb) 
{
	if (argc<idx+1){
		printf("Not enough arguments supplied, using %s as default\n", default_file);
	}
	else{
		*fb = fopen(argv[idx], "r");
		if (*fb != NULL) return 0;
		else{
			printf("Couldn't open file %s, using %s as default\n", argv[idx],default_file);
			fclose(*fb);
		}
	}
	*fb = fopen(default_file, "r");
	if (*fb == NULL){
		printf("Couldn't open %s\n",default_file);
		return -1;
	}
	return 0;
}

void freeData(Data * d, Cone * k){
	if(d) {
		if(d->b) coneOS_free(d->b);
		if(d->c) coneOS_free(d->c);
		if(d->Ax) coneOS_free(d->Ax);
		if(d->Ai) coneOS_free(d->Ai);
		if(d->Ap) coneOS_free(d->Ap);
		coneOS_free(d);
	}
	if(k) {
		if(k->q) coneOS_free(k->q);
		if(k->s) coneOS_free(k->s);
		coneOS_free(k);
	}
	d = NULL; k = NULL;
}

void freeSol(Sol *sol){
	if(sol) {
		if(sol->x) coneOS_free(sol->x);
		if(sol->y) coneOS_free(sol->y);
		//if(sol->status) coneOS_free(sol->status);
		coneOS_free(sol);
	}
	sol = NULL;
}


