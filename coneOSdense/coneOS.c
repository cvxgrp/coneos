/* coneos 1.0 */
#include "coneOS.h"
#include "normalize.h"

static int _lineLen_;
// constants and data structures
static const char* HEADER[] = {
	" Iter ", 
	"nm(u - u_t)",
	"nm(u - u_p)",
	"  kap/tau  ",
	"  time (s) ",
};

static const int HEADER_LEN = 5;

// to hold residual information
struct residuals {
	double resDual;
	double resPri;
	double tau;
	double kap;
};

static inline void updateDualVars(Data * d, Work * w);
static inline void projectCones(Data * d,Work * w,Cone * k);
static inline void sety(Data * d, Work * w, Sol * sol);
static inline void setx(Data * d, Work * w, Sol * sol);
static inline void sets(Data * d, Work * w, Sol * sol);
static inline int getSolution(Data * d, Work * w, Sol * sol, Info * info);
static inline void getInfo(Data * d, Work * w, Sol * sol, Info * info, struct residuals * r, int status);
static inline void printSummary(Data * d,Work * w,int i, struct residuals *r);
static inline void printHeader(Data * d, Work * w, Cone * k);
static inline void printFooter(Data * d, Info * info, int status); 
static inline void printSol(Data * d, Sol * sol, Info * info);
static inline void freeWork(Work * w);
static inline void projectLinSys(Data * d,Work * w);
static inline Work * initWork(Data * d, Cone * k);
static inline int converged(Data * d, Work * w, struct residuals * r);

/* coneOS returns the following integers:
   -2 failure
   -1 undetermined
   0 feasible, solved
   1 primal infeasible, dual unbounded
   2 primal unbounded, dual infeasible
 */

int coneOS(Data * d, Cone * k, Sol * sol, Info * info)
{
	if(d == NULL || k == NULL) {
		return -2;
	}
	tic();
	int i;
	struct residuals r = {-1, -1, -1};
	Work * w = initWork(d,k);
	if(d->VERBOSE) {
		printHeader(d, w, k);
	} /* coneOS: */
	for (i=0; i < d->MAX_ITERS; ++i){
		memcpy(w->u_prev, w->u, w->l*sizeof(double));
		
		projectLinSys(d,w);
		projectCones(d,w,k);
		updateDualVars(d,w);
	
		if (converged(d,w,&r)) break;

		if (i % 100 == 0){
			if (d->VERBOSE) printSummary(d,w,i,&r);
		}
	}
	if(d->VERBOSE) printSummary(d,w,i,&r);
	int status = getSolution(d,w,sol,info);

	if(d->NORMALIZE) unNormalize(d,w,sol,status);
	
	info->iter = i;
	getInfo(d,w,sol,info,&r,status);
	if(d->VERBOSE) printFooter(d, info, status);
	freeWork(w);
	return status;
}

static inline int converged(Data * d, Work * w, struct residuals * r){
	double tau = fabs((w->u[w->l-1]+w->u_t[w->l-1])/2);
	double kap = fabs(w->v[w->l-1]);
	r->resPri = calcNormDiff(w->u, w->u_t, w->l);
	r->resDual = calcNormDiff(w->u, w->u_prev, w->l);
    r->tau = tau;
    r->kap = kap;
	if (fmin(tau,kap)/fmax(tau,kap) < 1e-6 && fmax(r->resPri,r->resDual) < d->EPS_ABS*(tau+kap))
	{
		return 1;
	}
	return 0;
}

static inline void getInfo(Data * d, Work * w, Sol * sol, Info * info, struct residuals * r, int status){
	info->time = tocq();
 	
	double * dr = coneOS_calloc(d->n,sizeof(double));
	double * pr = coneOS_calloc(d->m,sizeof(double));
	double * x = sol->x, * y = sol->y, * s = sol->s;
	
	accumByA(d,x,pr); // pr = Ax
	addScaledArray(pr,s,d->m,1.0); // pr = Ax + s
	
	accumByAtrans(d,y,dr); // dr = A'y
	
    double cTx = innerProd(x,d->c,d->n);
	double bTy = innerProd(y,d->b,d->m);
    info->pobj = cTx;
	info->dobj = -bTy;

    if (status == 2) // dual infeasible
    {    
        info->dobj = NAN;
        info->gap = NAN;
        info->presid = calcNorm(pr,d->m);
        info->dresid = NAN;
    }
    else if (status == 1) // primal infeasible
    {
        info->pobj = NAN;
        info->gap = NAN;
        info->presid = NAN;
        info->dresid = calcNorm(dr,d->n);
    }
    else {
        addScaledArray(pr,d->b,d->m,-1.0); // pr = Ax + s - b
        addScaledArray(dr,d->c,d->n,1.0); // dr = A'y + c
        info->gap = info->pobj-info->dobj;
        info->presid = calcNorm(pr,d->m);
        info->dresid = calcNorm(dr,d->n);
    }
    coneOS_free(dr); coneOS_free(pr);
}



static inline Work * initWork(Data *d, Cone * k) {

	Work * w = coneOS_malloc(sizeof(Work));

	if(d->NORMALIZE) {
		normalize(d,w,k);
	}
	else {
		w->D = NULL;
		w->E = NULL;
		w->sc_c = 1.0;
		w->sc_b = 1.0;
		w->scale = 1.0;
	}

	w->l = d->n+d->m+1;
	w->u = coneOS_calloc(w->l,sizeof(double));
	w->u[w->l-1] = 1.0;
	w->v = coneOS_calloc(w->l,sizeof(double));
	//w->v[w->l-1] = 0.0;
	w->u_t = coneOS_calloc(w->l,sizeof(double));
	w->u_prev = coneOS_calloc(w->l,sizeof(double));
	w->h = coneOS_calloc((w->l-1),sizeof(double));
	memcpy(w->h,d->c,d->n*sizeof(double));
	memcpy(&(w->h[d->n]),d->b,d->m*sizeof(double));
	w->g = coneOS_calloc((w->l-1),sizeof(double));
	memcpy(w->g,w->h,(w->l-1)*sizeof(double));

	if (k->s){
		/* eigenvector decomp workspace */
		int i, nMax = 0;
		for (i=0; i < k->ssize; ++i){
			if (k->s[i] > nMax) nMax = k->s[i];
		}
		w->Xs = coneOS_calloc(nMax*nMax,sizeof(double));
		w->Z = coneOS_calloc(nMax*nMax,sizeof(double));
		w->e = coneOS_calloc(nMax,sizeof(double));
	} else {
		w->Xs = NULL;
		w->Z = NULL;
		w->e = NULL;
	}

	/* initialize the private data: */
	int status = privateInitWork(d, w);
	if (status < 0){
		coneOS_printf("privateInitWork failure: %i\n",status);
		exit(-1);
	} //else coneOS_printf("privateInitWork success: %i\n",status);
	d->CG_MAX_ITS = d->CG_MAX_ITS*100;
	d->CG_TOL = d->CG_TOL/100;
	solveLinSys(d,w,w->g, NULL); 
	d->CG_MAX_ITS = d->CG_MAX_ITS/100;
	d->CG_TOL = d->CG_TOL*100;
	scaleArray(&(w->g[d->n]),-1,d->m);
	w->gTh = innerProd(w->h, w->g, w->l-1); 
	return w;
}

static inline void projectLinSys(Data * d,Work * w){

	// ut = u + v
	memcpy(w->u_t,w->u,w->l*sizeof(double));
	addScaledArray(w->u_t,w->v,w->l,1.0);
	
	scaleArray(w->u_t,d->RHO_X,d->n);

	addScaledArray(w->u_t,w->h,w->l-1,-w->u_t[w->l-1]);
	addScaledArray(w->u_t, w->h, w->l-1, -innerProd(w->u_t,w->g,w->l-1)/(w->gTh+1));
	scaleArray(&(w->u_t[d->n]),-1,d->m);
	
	solveLinSys(d, w, w->u_t, w->u);
	
	w->u_t[w->l-1] += innerProd(w->u_t,w->h,w->l-1);
}

static inline void freeWork(Work * w){
	freePriv(w);
	if(w){
		if(w->Xs) coneOS_free(w->Xs);
		if(w->Z) coneOS_free(w->Z);
		if(w->e) coneOS_free(w->e);
		if(w->u) coneOS_free(w->u);
		if(w->v) coneOS_free(w->v);
		if(w->u_t) coneOS_free(w->u_t);
		if(w->u_prev) coneOS_free(w->u_prev);
		if(w->h) coneOS_free(w->h);
		if(w->g) coneOS_free(w->g);
		if(w->D) coneOS_free(w->D);
		if(w->E) coneOS_free(w->E);
		coneOS_free(w);
	}
}

static inline void printSol(Data * d, Sol * sol, Info * info){
	int i;
	coneOS_printf("%s\n",info->status); 
	if (sol->x != NULL){
		for ( i=0;i<d->n; ++i){
			coneOS_printf("x[%i] = %4f\n",i, sol->x[i]);
		}
	}
	if (sol->y != NULL){
		for ( i=0;i<d->m; ++i){
			coneOS_printf("y[%i] = %4f\n",i, sol->y[i]);
		}
	}
}

static inline void updateDualVars(Data * d, Work * w){
	int i;
	/*
	   for(i = 0; i < d->n; ++i) { 
	   w->v[i] += w->u[i] - w->u_t[i]; 
	   }
	 */
	//for(i = 0; i < w->l; ++i) { 
	if (fabs(d->ALPH - 1.0) < 1e-9) {
		// this is over-step parameter:
		//double sig = (1+sqrt(5))/2;
		double sig = 1.0;
		for(i = d->n; i < w->l; ++i) { 
			w->v[i] += sig*(w->u[i] - w->u_t[i]);
		}
	}
	else {
		// this does not relax 'x' variable
		for(i = d->n; i < w->l; ++i) { 
			w->v[i] += (w->u[i] - d->ALPH*w->u_t[i] - (1.0 - d->ALPH)*w->u_prev[i]); 
		}
	}
}

static inline void projectCones(Data *d,Work * w,Cone * k){
	int i;
	// this does not relax 'x' variable
	for(i = 0; i < d->n; ++i) { 
		w->u[i] = w->u_t[i] - w->v[i];
	}
	//for(i = 0; i < w->l; ++i){
	for(i = d->n; i < w->l; ++i){
		w->u[i] = d->ALPH*w->u_t[i] + (1-d->ALPH)*w->u_prev[i] - w->v[i];
	}
	/* u = [x;y;tau] */
	projCone(&(w->u[d->n]),k,w);
	if (w->u[w->l-1]<0.0) w->u[w->l-1] = 0.0;
}

static inline int getSolution(Data * d,Work * w,Sol * sol, Info * info){
    double tau = (w->u[w->l-1]+w->u_t[w->l-1])/2;
    //double tau = w->u[w->l-1];
    double kap = fabs(w->v[w->l-1]);
    setx(d,w,sol);
    sety(d,w,sol);
    sets(d,w,sol);
    if (tau > d->UNDET_TOL && tau > kap){
        memcpy(info->status,"Solved", 7);
        scaleArray(sol->x,1.0/tau,d->n);
        scaleArray(sol->y,1.0/tau,d->m);
        scaleArray(sol->s,1.0/tau,d->m);
        return 0;
    }   
    else{ 
        if (calcNorm(w->u,w->l)<d->UNDET_TOL*sqrt(w->l)){
            memcpy(info->status, "Indeterminate", 15);
            scaleArray(sol->x,NAN,d->n);
            scaleArray(sol->y,NAN,d->m);
            scaleArray(sol->s,NAN,d->m);
            return -1;
        }   
        else {
            double ip_y = innerProd(d->b,sol->y,d->m);  
            double ip_x = innerProd(d->c,sol->x,d->n);  
            if (ip_y < ip_x){
                memcpy(info->status,"Infeasible", 12);
                //scaleArray(sol->y,-1/ip_y,d->m);
                scaleArray(sol->x,NAN,d->n);
                scaleArray(sol->s,NAN,d->m);
                return 1;
            }   
            else{
                memcpy(info->status,"Unbounded", 11);
                //scaleArray(sol->x,-1/ip_x,d->n);
                scaleArray(sol->y,NAN,d->m);
                //scaleArray(sol->s,0.0,d->m);
                return 2;
            }
        }
    }
}

static inline void sety(Data * d,Work * w, Sol * sol){
	sol->y = coneOS_malloc(sizeof(double)*d->m);
	int i;
	for(i = 0; i < d->m; ++i) {
		/* not taking average has desirable properties, like returning sparse
		   solutions if l1 penalty used (for example), but will have worse accuracy
		   and worse dual equality constraint residual */
		   //sol->y[i] = 0.5 * (w->u[i + d->n]+w->u_t[i + d->n]);
		   sol->y[i] = w->u[i + d->n];
	}
}

static inline void sets(Data * d,Work * w, Sol * sol){
    sol->s = coneOS_malloc(sizeof(double)*d->m);
    memcpy(sol->s, &(w->v[d->n]), d->m*sizeof(double));
}

static inline void setx(Data * d,Work * w, Sol * sol){
	sol->x = coneOS_malloc(sizeof(double)*d->n);
	memcpy(sol->x, w->u, d->n*sizeof(double));
}

static inline void printSummary(Data * d,Work * w,int i, struct residuals *r){
	coneOS_printf("%*i|", (int)strlen(HEADER[0]), i);
	coneOS_printf("%*.2e ", (int)strlen(HEADER[1])-1, r->resPri/(r->tau+r->kap));
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[2])-1, r->resDual/(r->tau+r->kap));
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[3])-1, r->kap/r->tau);
   	coneOS_printf(" %*.2e ", (int)strlen(HEADER[4])-1, tocq()/1e3);
    coneOS_printf("\n");
#ifdef MATLAB_MEX_FILE
	mexEvalString("drawnow;");
#endif
}

static inline void printHeader(Data * d, Work * w, Cone * k) {
	int i;  
	_lineLen_ = -1;
	for(i = 0; i < HEADER_LEN; ++i) {
		_lineLen_ += strlen(HEADER[i]) + 1;
	}
	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
	coneOS_printf("\nconeOS 1.0: %s method\n", w->method);
   	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
	coneOS_printf("\nvariables n = %i, constraints m = %i\n", d->n, d->m);
    coneOS_printf("EPS = %.2e\n", d->EPS_ABS);
    coneOS_printf("cones:\tfree/zero %i\n\tlinear %i\n\tsecond-order %i\n\tsemi-definite %i\n", k->f, k->l, k->qsize, k->ssize);
    
    for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
	coneOS_printf("\n");
	for(i = 0; i < HEADER_LEN - 1; ++i) {
		coneOS_printf("%s|", HEADER[i]);
	}
	coneOS_printf("%s\n", HEADER[HEADER_LEN-1]);
	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("=");
	}
	coneOS_printf("\n");
}

static inline void printFooter(Data * d, Info * info, int status) {
    double b_inf = calcNormInf(d->b, d->m);
    double c_inf = calcNormInf(d->c, d->n);
    double gap_rel = 1 + fabs(info->pobj) + fabs(info->dobj);
	int i;  
	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
	coneOS_printf("\nStatus: %s\n",info->status);
	if (info->iter == d->MAX_ITERS)
		coneOS_printf("Hit MAX_ITERS, solution may be inaccurate\n"); 
	//coneOS_printf("Primal objective value: %4f\n",info->pobj);
	//coneOS_printf("Dual objective value: %4f\n",info->dobj);
	//coneOS_printf("Duality gap: %e\n", info->gap);
	coneOS_printf("Time taken: %.4f seconds\n",info->time/1e3);

	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
    coneOS_printf("\n");

    if (status == 1) {
        coneOS_printf("Certificate of primal infeasibility:\n");
        coneOS_printf("|A'y|_2/(1+|c|_inf) = %2e\n", info->dresid/(1+c_inf));
        coneOS_printf("dist(y,K*) = 0\n");
        coneOS_printf("b'y = %.4f\n", -info->dobj);
    } 
    else if (status == 2) {
        coneOS_printf("Certificate of dual infeasibility:\n");
        coneOS_printf("|Ax+s|_2/(1+|b|_inf) = %2e\n", info->presid/(1+b_inf));
        coneOS_printf("dist(s,K) = 0\n");
        coneOS_printf("c'x = %.4f\n", info->pobj);
    }
    else {
        coneOS_printf("DIMACS error metrics:\n");
        coneOS_printf("|Ax+s-b|_2/(1+|b|_inf) = %2e\n|A'y+c|_2/(1+|c|_inf) = %2e\n",info->presid/(1+b_inf), info->dresid/(1+c_inf));
        coneOS_printf("|c'x+b'y|/(1+|c'x|+|b'y|) = %2e\n", fabs(info->gap/gap_rel)); 
        coneOS_printf("dist(s,K) = 0, dist(y,K*) = 0, s'y = 0\n");
        for(i = 0; i < _lineLen_; ++i) {
            coneOS_printf("-");
        }
        coneOS_printf("\n");
        coneOS_printf("c'x = %.4f, -b'y = %.4f\n",info->pobj, info->dobj);
    }
    for(i = 0; i < _lineLen_; ++i) {
        coneOS_printf("=");
    }
    coneOS_printf("\n");
}
