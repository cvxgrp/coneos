/* coneos 1.0 */
#include "coneOS.h"

static int _lineLen_;
// constants and data structures
static const char* HEADER[] = {
	" Iter ", 
	"nm(u - u_t)",
	"nm(u - u_p)",
	"  pri res  ",
	"  dual res ",
	"  pri obj  ",
	"  dual obj ",
	"    gap    ",
	"  kap/tau  ",
};

static const int HEADER_LEN = 9;

// to hold residual information
struct residuals {
	double resDual;
	double resPri;
	double pres;
	double dres;
	double pobj;
	double dobj;
	double dgap;
	double tau;
	double kap;
};

static inline void updateDualVars(Data * d, Work * w);
static inline void projectCones(Data * d,Work * w,Cone * k);
static inline void sety(Data * d, Work * w, Sol * sol);
static inline void setx(Data * d, Work * w, Sol * sol);
static inline int getSolution(Data * d, Work * w, Sol * sol, Info * info);
static inline void getInfo(Data * d, Work * w, Sol * sol, Info * info, struct residuals * r, int status);
static inline void printSummary(Data * d,Work * w,int i, struct residuals *r);
static inline void printHeader(Data * d, Work * w);
static inline void printFooter(Data * d, Info * info); 
static inline void printSol(Data * d, Sol * sol, Info * info);
static inline void freeWork(Work * w);
static inline void projectLinSys(Data * d,Work * w);
static inline Work * initWork(Data * d, Cone * k);
static inline void unNormalize(Data *d, Work * w, Sol * sol);
static inline int converged(Data * d, Work * w, struct residuals * r);

/* coneOS returns the following integers:
   -2 failure
   -1 undetermined
   0 feasible, solved
   1 infeasible
   2 unbounded
 */

int coneOS(Data * d, Cone * k, Sol * sol, Info * info)
{
	if(d == NULL || k == NULL) {
		return -2;
	}
	tic();
	int i;
	struct residuals r = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
	Work * w = initWork(d,k);
	if(d->VERBOSE) {
		printHeader(d, w);
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

	if(d->NORMALIZE) unNormalize(d,w,sol);
	
	info->iter = i;
	getInfo(d,w,sol,info,&r,status);
	if(d->VERBOSE) printFooter(d, info);
	freeWork(w);
	return status;
}

static inline int converged(Data * d, Work * w, struct residuals * r){
	double tau = fabs((w->u[w->l-1]+w->u_t[w->l-1])/2);
	double kap = fabs(w->v[w->l-1]);
	r->resPri = calcNormDiff(w->u, w->u_t, w->l);
	r->resDual = calcNormDiff(w->u, w->u_prev, w->l);
	
	if (fmin(tau,kap)/fmax(tau,kap) < 1e-6 && fmax(r->resPri,r->resDual) < d->EPS_ABS*(tau+kap))
	{
		return 1;
	}
	return 0;
}

static inline void getInfo(Data * d, Work * w, Sol * sol, Info * info, struct residuals * r,int status){
	info->time = tocq();
	info->presid = r->pres;
	info->dresid = r->dres;
	info->pobj = r->pobj;
	info->dobj = r->dobj;
	if (status == 2)
		info->pobj = NAN;
	else if(status == 1)
		info->dobj = NAN;
	info->gap = info->pobj-info->dobj;
}

static inline void unNormalize(Data *d, Work * w, Sol * sol){
	int i, j;
	double * D = w->D;
	double * E = w->E;
	for (i = 0; i < d->n; ++i){
		sol->x[i] /= (w->E[i] * w->sc_b);
	}
	for (i = 0; i < d->m; ++i){
		sol->y[i] /= (w->D[i] * w->sc_c);
	}
	for (i = 0; i < d->n; ++i){
		d->c[i] *= E[i]/(w->sc_c * w->scale);
	}
	for (i = 0; i < d->m; ++i){
		d->b[i] *= D[i]/(w->sc_b * w->scale);
	}
    for(i = 0; i < d->n; ++i){
        for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
            d->Ax[j] *= D[d->Ai[j]];
        }   
    }   
    for (i = 0; i < d->n; ++i){
        scaleArray(&(d->Ax[d->Ap[i]]), E[i]/w->scale, d->Ap[i+1] - d->Ap[i]);
    }   
}

static inline void normalize(Data * d, Work * w, Cone * k){
    
	double * D = coneOS_calloc(d->m, sizeof(double));
	double * E = coneOS_malloc(d->n*sizeof(double));
	
	int i, j, count;
	double wrk;
	
	// heuristic rescaling, seems to do well with a scaling of about 4
	w->scale = 4.0;

	// calculate row norms
	for(i = 0; i < d->n; ++i){
		for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
			wrk = d->Ax[j];
			D[d->Ai[j]] += wrk*wrk;
		}
	}
	for (i=0; i < d->m; ++i){
		D[i] = fmax(sqrt(D[i]),1e-6); // just the norms
	}
	// mean of norms of rows across each cone	
    count = k->l+k->f;
	for(i = 0; i < k->qsize; ++i)
    {
		wrk = 0;
		for (j = count; j < count + k->q[i]; ++j){
        	wrk += D[j];
		}
		wrk /= k->q[i];
		for (j = count; j < count + k->q[i]; ++j){
        	D[j] = wrk;
		}
		count += k->q[i];
    }
    for (i=0; i < k->ssize; ++i)
	{
 		wrk = 0;
		for (j = count; j < count + k->s[i]; ++j){
        	wrk += D[j];
		}
		wrk /= k->s[i];
		for (j = count; j < count + k->s[i]; ++j){
        	D[j] = wrk;
		}
		count += (k->s[i])*(k->s[i]);
    }
	// scale the rows with D
	for(i = 0; i < d->n; ++i){
		for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
			d->Ax[j] /= D[d->Ai[j]];
		}
	}
	// calculate and scale by col norms, E
	for (i = 0; i < d->n; ++i){
		E[i] = fmax(calcNorm(&(d->Ax[d->Ap[i]]),d->Ap[i+1] - d->Ap[i]),1e-6);
		scaleArray(&(d->Ax[d->Ap[i]]), w->scale/E[i], d->Ap[i+1] - d->Ap[i]);
	}
	// scale b
	for (i = 0; i < d->m; ++i){
		d->b[i] /= D[i];
	}
	w->sc_b = 1/fmax(calcNorm(d->b,d->m),1e-6);
	scaleArray(d->b, w->scale * w->sc_b, d->m);
	// scale c
	for (i = 0; i < d->n; ++i){
		d->c[i] /= E[i];
	}
	double meanNormRowA = 0.0;
	double *nms = coneOS_calloc(d->m,sizeof(double));
	for(i = 0; i < d->n; ++i){
		for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
			wrk = d->Ax[j];
			nms[d->Ai[j]] += wrk*wrk;
		}
	}
	for (i=0; i < d->m; ++i){
		meanNormRowA += sqrt(nms[i])/d->m;
	}
	w->sc_c = meanNormRowA/fmax(calcNorm(d->c,d->n),1e-6);
	scaleArray(d->c, w->scale * w->sc_c, d->n);

	w->D = D;
	w->E = E;

	coneOS_free(nms);
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
	w->v[w->l-1] = 0.0;
	w->u_t = coneOS_malloc(w->l*sizeof(double));
	w->u_prev = coneOS_malloc(w->l*sizeof(double));
	w->h = coneOS_malloc((w->l-1)*sizeof(double));
	memcpy(w->h,d->c,d->n*sizeof(double));
	memcpy(&(w->h[d->n]),d->b,d->m*sizeof(double));
	w->g = coneOS_malloc((w->l-1)*sizeof(double));
	memcpy(w->g,w->h,(w->l-1)*sizeof(double));

	if (k->s){
		/* eigenvector decomp workspace */
		int i, nMax = 0;
		for (i=0; i < k->ssize; ++i){
			if (k->s[i] > nMax) nMax = k->s[i];
		}
		w->Xs = coneOS_malloc(nMax*nMax*sizeof(double));
		w->Z = coneOS_malloc(nMax*nMax*sizeof(double));
		w->e = coneOS_malloc(nMax*sizeof(double));
	} else {
		w->Xs = NULL;
		w->Z = NULL;
		w->e = NULL;
	}

	/* initialize the private data: */
	privateInitWork(d, w);
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
	if (tau > d->UNDET_TOL && tau > kap){
		memcpy(info->status,"Solved", 7);
		scaleArray(sol->x,1.0/tau,d->n);
		scaleArray(sol->y,1.0/tau,d->m);
		return 0;
	}   
	else{ 
		if (calcNorm(w->u,w->l)<d->UNDET_TOL*sqrt(w->l)){
			memcpy(info->status, "Indeterminate", 15);
			scaleArray(sol->x,NAN,d->n);
			scaleArray(sol->y,NAN,d->m);
			return -1;
		}   
		else {
			double ip_y = innerProd(d->b,sol->y,d->m);  
			double ip_x = innerProd(d->c,sol->x,d->n);  
			if (ip_y < ip_x){
				memcpy(info->status,"Infeasible", 12);
				scaleArray(sol->y,-1/ip_y,d->m);
				scaleArray(sol->x,NAN,d->n);
				return 1;
			}   
			else{
				memcpy(info->status,"Unbounded", 11);
				scaleArray(sol->x,-1/ip_x,d->n);
				scaleArray(sol->y,NAN,d->m);
				return 2;
			}   
		}   
	}   
}

static inline void sety(Data * d,Work * w, Sol * sol){
	sol->y = coneOS_malloc(sizeof(double)*d->m);
	int i;
	for(i = 0; i < d->m; ++i) {
		/* XXX: (see also converged function) 
		   not taking average has desirable properties, like returning sparse
		   solutions if l1 penalty used (for example), but will have worse accuracy
		   and worse dual equality constraint residual */
		sol->y[i] = 0.5 * (w->u[i + d->n]+w->u_t[i + d->n]);
		//sol->y[i] = w->u[i + d->n];
	}
}

static inline void setx(Data * d,Work * w, Sol * sol){
	sol->x = coneOS_malloc(sizeof(double)*d->n);
	memcpy(sol->x, w->u, d->n*sizeof(double));
}

static inline void printSummary(Data * d,Work * w,int i, struct residuals *r){
	// coneOS_printf("Iteration %i, primal residual %4f, primal tolerance %4f\n",i,err,EPS_PRI);
	double tau = fabs(w->u[w->l-1]+w->u_t[w->l-1])/2;
	//double tau = w->u[w->l-1];
	double kap = fabs(w->v[w->l-1]);
 	
	double * dr = coneOS_calloc(d->n,sizeof(double));
	double * pr = coneOS_calloc(d->m,sizeof(double));

	/* see discussion in sety function
	double y[d->m], * x = w->u;
	memcpy(y,&(w->u[d->n]),d->m*sizeof(double));
	addScaledArray(y,&(w->u_t[d->n]),d->m,1);
	scaleArray(y,0.5,d->m);
	*/
	double * x = w->u, * y = &(w->u[d->n]);
	
	// slower primal resid calculation (mult by A)
	/*		
	accumByA(d,x,pr);
	addScaledArray(pr,&(w->v[d->n]),d->m,1.0);
	addScaledArray(pr,d->b,d->m,-tau);
	*/
	
	// better primal resid calculation
	double * y_prev = &(w->u_prev[d->n]);
	double * y_t = &(w->u_t[d->n]);
	double dt = w->u[w->l-1] - (d->ALPH*w->u_t[w->l-1] + (1-d->ALPH)*(w->u_prev[w->l-1]));
	int j;
	for (j = 0; j < d->m; ++j){
		pr[j] = d->b[j]*dt - y[j] + (2-d->ALPH)*y_prev[j] - (1-d->ALPH)*y_t[j];	
	}
	
	r->pres = calcNorm(pr,d->m)/tau;

	accumByAtrans(d,y,dr);
	addScaledArray(dr,d->c,d->n,tau);
	r->dres = calcNorm(dr,d->n)/tau;
	
	double cTx = innerProd(x,d->c,d->n);
	double bTy = innerProd(y,d->b,d->m);

	r->dgap = fabs(cTx + bTy);
	r->pobj = cTx;
	r->dobj = -bTy;
	r->tau = tau;
	r->kap = kap;
	coneOS_free(dr); coneOS_free(pr);

	coneOS_printf("%*i|", (int)strlen(HEADER[0]), i);
	coneOS_printf("%*.2e ", (int)strlen(HEADER[1])-1, r->resPri);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[2])-1, r->resDual);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[3])-1, r->pres); 
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[4])-1, r->dres);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[5])-1, r->pobj);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[6])-1, r->dobj);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[7])-1, r->dgap);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[8])-1, r->kap/r->tau);
coneOS_printf("\n");
#ifdef MATLAB_MEX_FILE
	mexEvalString("drawnow;");
#endif
}

static inline void printHeader(Data * d, Work * w) {
	int i;  
	_lineLen_ = -1;
	for(i = 0; i < HEADER_LEN; ++i) {
		_lineLen_ += strlen(HEADER[i]) + 1;
	}

	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
	coneOS_printf("\nconeOS 1.0: %s method, A matrix density: %4f\n",w->method,((double)d->Anz/d->n)/d->m);
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

static inline void printFooter(Data * d, Info * info) {
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
	coneOS_printf("Time taken: %4f seconds\n",info->time/1e3);

	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("=");
	}
	coneOS_printf("\n");
}
