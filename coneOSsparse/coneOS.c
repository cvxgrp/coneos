/* coneos 1.0 */
#include "coneOS.h"
#include "normalize.h"

static int _lineLen_;
// constants and data structures
static const char* HEADER[] = {
	" Iter ", 
	" pri res ",
	" dua res ",
    " rel gap ",
	" pri obj ",
    " dua obj ",
    "  kappa  ",
	" time (s)",
};

static const int HEADER_LEN = 8;

static inline void updateDualVars(Data * d, Work * w);
static inline void projectCones(Data * d,Work * w,Cone * k);
static inline void sety(Data * d, Work * w, Sol * sol);
static inline void setx(Data * d, Work * w, Sol * sol);
static inline void sets(Data * d, Work * w, Sol * sol);
static inline void setSolution(Data * d, Work * w, Sol * sol, Info * info);
static inline void getInfo(Data * d, Work * w, Sol * sol, Info * info, struct residuals * r);
static inline void printSummary(Data * d,Work * w,int i, struct residuals *r);
static inline void printHeader(Data * d, Work * w, Cone * k);
static inline void printFooter(Data * d, Info * info, Work * w); 
static inline void freeWork(Work * w);
static inline void projectLinSys(Data * d,Work * w);
static inline Work * initWork(Data * d, Cone * k);
static inline int converged(Data * d, Work * w, struct residuals * r, int iter);
static inline int exactConverged(Data * d, Work * w, struct residuals * r, int iter);

#define PRINT_INTERVAL 100
#define CONVERGED_INTERVAL 20

/* coneOS returns one of the following integers: */
/* (zero should never be returned) */

#define FAILURE -4
#define INDETERMINATE -3
#define INFEASIBLE -2 // primal infeasible, dual unbounded
#define UNBOUNDED -1 // primal unbounded, dual infeasible
#define SOLVED 1

int coneOS(Data * d, Cone * k, Sol * sol, Info * info)
{
	if(d == NULL || k == NULL) {
		return FAILURE;
	}
    tic();
	info->stint = 0; // not yet converged
    int i;
	struct residuals r = {-1, -1, -1, -1, -1, -1, -1};
    Work * w = initWork(d,k);
	if(d->VERBOSE) {
		printHeader(d, w, k);
	} /* coneOS: */
	for (i=0; i < d->MAX_ITERS; ++i){
		memcpy(w->u_prev, w->u, w->l*sizeof(double));
		
		projectLinSys(d,w);
		projectCones(d,w,k);
		updateDualVars(d,w);
	    
        info->stint = converged(d,w,&r,i);

		if (info->stint != 0) break;

		if (i % PRINT_INTERVAL == 0){
			if (d->VERBOSE) printSummary(d,w,i,&r);
		}
	}
	if(d->VERBOSE) printSummary(d,w,i,&r);
	setSolution(d,w,sol,info);

	if(d->NORMALIZE) unNormalize(d,w,sol);

    info->iter = i;
	
    getInfo(d,w,sol,info,&r);
	if(d->VERBOSE) printFooter(d, info, w);
	freeWork(w);
	return info->stint;
}

static inline int converged(Data * d, Work * w, struct residuals * r, int iter){
    /* approximate convergence check:
    double tau = fabs(w->u[w->l-1]); // abs to prevent negative stopping tol
    double kap = fabs(w->v[w->l-1]);
    r->resPri = calcNormDiff(w->u, w->u_t, w->l);
    r->resDual = calcNormDiff(w->u, w->u_prev, w->l);
    r->tau = tau;
    r->kap = kap;
    if (fmin(tau,kap)/fmax(tau,kap) < 1e-6 && fmax(r->resPri, r->resDual) < d->EPS_ABS*(tau+kap)){
        return 1;
    }
    */
    if (iter % CONVERGED_INTERVAL == 0) {
        return exactConverged(d,w,r,iter);
    }
    return 0;
}

static inline int exactConverged(Data * d, Work * w, struct residuals * r, int iter){
    double * pr = coneOS_calloc(d->m,sizeof(double));
    double * dr = coneOS_calloc(d->n,sizeof(double));
    double * Axs = coneOS_calloc(d->m,sizeof(double));
    double * ATy = coneOS_calloc(d->n,sizeof(double));

    double tau = fabs(w->u[w->l-1]);
    double kap = fabs(w->v[w->l-1]);
    double * x = w->u, * y = &(w->u[d->n]);
    double * D = w->D, * E = w->E;
    int i;

    /*
    // requires mult by A:
    double * s = &(w->v[d->n]);
    accumByA(d,x,Axs); // Axs = Ax
    addScaledArray(Axs,s,d->m,1.0); // Axs = Ax + s
    memcpy(pr, Axs, d->m * sizeof(double)); // pr = Ax + s
    addScaledArray(pr,d->b,d->m,-tau); // pr = Ax + s - b * tau
    */

    // does not require mult by A:
    memcpy(pr,&(w->u[d->n]),d->m * sizeof(double));
    addScaledArray(pr,&(w->u_prev[d->n]),d->m,d->ALPH-2);
    addScaledArray(pr,&(w->u_t[d->n]),d->m,1-d->ALPH);
    addScaledArray(pr,d->b, d->m, w->u_t[w->l-1] - tau) ; // pr = Ax + s - b * tau
    memcpy(Axs, pr, d->m * sizeof(double));
    addScaledArray(Axs, d->b, d->m, tau); // Axs = Ax + s

    double cTx = innerProd(x,d->c,d->n);

    if (d->NORMALIZE) {
        kap /= (w->scale * w->sc_c * w->sc_b);
        for (i = 0; i < d->m; ++i) {
            pr[i] *= D[i]/(w->sc_b * w->scale);
            Axs[i] *= D[i]/(w->sc_b * w->scale);
        } 
        cTx /= (w->scale * w->sc_c * w->sc_b);
    }
    r->tau = tau;
    r->kap = kap;

    double nmAxs = calcNorm(Axs,d->m);
    r->resPri = cTx < 0 ? w->nm_c * nmAxs / -cTx : NAN;
    //coneOS_printf("unbounded cert: %4e\n", w->nm_c * nmAxs / (1+w->nm_b) / -cTx);
    if (r->resPri < d->EPS_ABS) {
        return UNBOUNDED;
    }

    accumByAtrans(d,y,ATy); // ATy = A'y
    memcpy(dr, ATy, d->n * sizeof(double));
    addScaledArray(dr,d->c,d->n,tau); // dr = A'y + c * tau    

    double bTy = innerProd(y,d->b,d->m);

    if (d->NORMALIZE) {
        for (i = 0; i < d->n; ++i) {
            dr[i] *= E[i]/(w->sc_c * w->scale);
            ATy[i] *= E[i]/(w->sc_c * w->scale);
        }
        bTy /= (w->scale * w->sc_c * w->sc_b);
    }

    double nmATy = calcNorm(ATy,d->n);
    r->resDual = bTy < 0 ? w->nm_b * nmATy / -bTy : NAN;
    //coneOS_printf("infeas cert: %4e\n", w->nm_b * nmATy / (1+w->nm_c) /  - bTy );
    if (r->resDual < d->EPS_ABS) {
        return INFEASIBLE;
    }
    r->relGap = NAN;

    int status = 0;
    if (tau > kap) {
        double rpri = calcNorm(pr,d->m) / (1+w->nm_b) / tau;
        double rdua = calcNorm(dr,d->n) / (1+w->nm_c) / tau;
        double gap = fabs(cTx + bTy) / (tau + fabs(cTx) + fabs(bTy));

        r->resPri = rpri;
        r->resDual = rdua;
        r->relGap = gap;
        r->cTx = cTx / tau;
        r->bTy = bTy / tau;
        // coneOS_printf("primal resid: %4e, dual resid %4e, pobj %4e, dobj %4e, gap %4e\n", rpri,rdua,cTx,-bTy,gap);
        // coneOS_printf("primal resid: %4e, dual resid %4e, gap %4e\n",rpri,rdua,gap);
        if (fmax(fmax(rpri,rdua),gap) < d->EPS_ABS) {
            status = SOLVED;
        }
    } else {
        r->cTx = NAN;
        r->bTy = NAN;
    }
    coneOS_free(dr); coneOS_free(pr); coneOS_free(Axs); coneOS_free(ATy);
    return status;
}

static inline void getInfo(Data * d, Work * w, Sol * sol, Info * info, struct residuals * r){
    double * x = sol->x, * y = sol->y, * s = sol->s;
    
    double * dr = coneOS_calloc(d->n,sizeof(double));
    double * pr = coneOS_calloc(d->m,sizeof(double));

    accumByA(d,x,pr); // pr = Ax
    addScaledArray(pr,s,d->m,1.0); // pr = Ax + s

    accumByAtrans(d,y,dr); // dr = A'y

    double cTx = innerProd(x,d->c,d->n);
    double bTy = innerProd(y,d->b,d->m);
    info->pobj = cTx;
    info->dobj = -bTy;
    if (info->stint == SOLVED){
        addScaledArray(pr,d->b,d->m,-1.0); // pr = Ax + s - b
        addScaledArray(dr,d->c,d->n,1.0); // dr = A'y + c
        info->relGap = fabs(cTx + bTy) / (1 + fabs(cTx) + fabs(bTy));
        info->resPri = calcNorm(pr,d->m) / (1 + w->nm_b);
        info->resDual = calcNorm(dr,d->n) / (1+ w->nm_c);
    } else {
        if (info->stint == UNBOUNDED) {    
            info->dobj = NAN;
            info->relGap = NAN;
            info->resPri = w->nm_c * calcNorm(pr,d->m) / -cTx ;
            info->resDual = NAN;
            scaleArray(x,-1/cTx,d->n);
            scaleArray(s,-1/cTx,d->m);
            info->pobj = -1;
        }
        else {
            info->pobj = NAN;
            info->relGap = NAN;
            info->resPri = NAN;
            info->resDual = w->nm_b * calcNorm(dr,d->n) / -bTy ;
            scaleArray(y,-1/bTy,d->m);
            info->dobj = -1;
        }
    }
    info->time = tocq();
    coneOS_free(dr); coneOS_free(pr);
}

static inline Work * initWork(Data *d, Cone * k) {

	Work * w = coneOS_malloc(sizeof(Work));
    
    w->nm_b = calcNorm(d->b, d->m);
    w->nm_c = calcNorm(d->c, d->n);

    //w->nm_b = calcNormInf(d->b, d->m);
    //w->nm_c = calcNormInf(d->c, d->n);
    //w->nm_Q = calcNormFroQ(d);

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
	w->u[w->l-1] = sqrt(w->l);
	w->v = coneOS_calloc(w->l,sizeof(double));
	w->v[w->l-1] = sqrt(w->l);
	w->u_t = coneOS_calloc(w->l,sizeof(double));
	w->u_prev = coneOS_calloc(w->l,sizeof(double));
	w->h = coneOS_calloc((w->l-1),sizeof(double));
	memcpy(w->h,d->c,d->n*sizeof(double));
	memcpy(&(w->h[d->n]),d->b,d->m*sizeof(double));
	w->g = coneOS_calloc((w->l-1),sizeof(double));
	memcpy(w->g,w->h,(w->l-1)*sizeof(double));
	/* initialize the private data: */
	int status = privateInitWork(d, w);
	if (status < 0){
		coneOS_printf("privateInitWork failure: %i\n",status);
		exit(-1);
	} 
    //else coneOS_printf("privateInitWork success: %i\n",status);
	
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
        if(w->method) coneOS_free(w->method);
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

void printSol(Data * d, Sol * sol, Info * info){
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

static inline int solved(Data * d, Sol * sol, Info * info, double tau){
    memcpy(info->status,"Solved", 7);
    scaleArray(sol->x,1.0/tau,d->n);
    scaleArray(sol->y,1.0/tau,d->m);
    scaleArray(sol->s,1.0/tau,d->m);
    return SOLVED;
}

static inline int indeterminate(Data * d, Sol * sol, Info * info){
    memcpy(info->status, "Indeterminate", 15);
    scaleArray(sol->x,NAN,d->n);
    scaleArray(sol->y,NAN,d->m);
    scaleArray(sol->s,NAN,d->m);
    return INDETERMINATE;
}

static inline int infeasible(Data * d, Sol * sol, Info * info){
    memcpy(info->status,"Infeasible", 12);
    //scaleArray(sol->y,-1/ip_y,d->m);
    scaleArray(sol->x,NAN,d->n);
    scaleArray(sol->s,NAN,d->m);
    return INFEASIBLE;
}

static inline int unbounded(Data * d, Sol * sol, Info * info){
    memcpy(info->status,"Unbounded", 11);
    //scaleArray(sol->x,-1/ip_x,d->n);
    scaleArray(sol->y,NAN,d->m);
    return UNBOUNDED;
}

static inline void setSolution(Data * d,Work * w,Sol * sol, Info * info){
    setx(d,w,sol);
    sety(d,w,sol);
    sets(d,w,sol);
    if (info->stint == 0 || info->stint == SOLVED){
        double tau = w->u[w->l-1];
        double kap = fabs(w->v[w->l-1]);
        if (tau > d->UNDET_TOL && tau > kap){
            info->stint = solved(d,sol,info,tau);
        }   
        else{ 
            if (calcNorm(w->u,w->l)<d->UNDET_TOL*sqrt(w->l)){
                info->stint = indeterminate(d,sol,info);
            }   
            else {
                double bTy = innerProd(d->b,sol->y,d->m);  
                double cTx = innerProd(d->c,sol->x,d->n);  
                if (bTy < cTx){
                    info->stint = infeasible(d,sol,info);
                }   
                else{
                    info->stint = unbounded(d,sol,info);
                }
            }
        }
    } else if (info->stint == INFEASIBLE) {
        info->stint = infeasible(d,sol,info);
    } else {
        info->stint = unbounded(d,sol,info);
    }
}

static inline void sety(Data * d,Work * w, Sol * sol){
	sol->y = coneOS_malloc(sizeof(double)*d->m);
	memcpy(sol->y, &(w->u[d->n]), d->m*sizeof(double));
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
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[1])-1, r->resPri);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[2])-1, r->resDual);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[3])-1, r->relGap);
    if (r->cTx < 0) {
	    coneOS_printf("%*.2e ", (int)strlen(HEADER[4])-1, r->cTx);
    } else {
        coneOS_printf(" %*.2e ", (int)strlen(HEADER[4])-1, r->cTx);
    }
    if (r->bTy >= 0) {
	    coneOS_printf("%*.2e ", (int)strlen(HEADER[5])-1, -r->bTy);
    } else {
	    coneOS_printf(" %*.2e ", (int)strlen(HEADER[5])-1, -r->bTy);
    }
    coneOS_printf(" %*.2e ", (int)strlen(HEADER[6])-1, r->kap);
	coneOS_printf(" %*.2e ", (int)strlen(HEADER[7])-1, tocq()/1e3);
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
	coneOS_printf("\nconeOS 1.0: %s\n", w->method);
   	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
    coneOS_printf("\nEPS = %.2e, ALPHA = %.2f, MAX_ITERS = %i, NORMALIZE = %i\n", d->EPS_ABS, d->ALPH, d->MAX_ITERS, d->NORMALIZE);
	coneOS_printf("variables n = %i, constraints m = %i, non-zeros in A = %i\n", d->n, d->m, d->Anz);
    
    int socVars = 0;
    for (int i=0;i<k->qsize;i++){
        socVars += k->q[i];
    }
    int sdVars = 0;
    for (int i=0;i<k->ssize;i++){
        sdVars += k->s[i]*k->s[i];
    }

    coneOS_printf("cones:\tzero/free vars: %i\n\tlinear vars: %i\n\tsoc vars: %i, soc blks: %i\n\tsd vars: %i, sd blks: %i\n\texp vars: %i\n\tdual exp vars: %i\n", k->f, k->l, socVars, k->qsize, sdVars,k->ssize, k->ep*3, k->ed*3);
    
 
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

static inline void printFooter(Data * d, Info * info, Work * w) {
	int i;
	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
	coneOS_printf("\nStatus: %s\n",info->status);
	if (info->iter == d->MAX_ITERS) {
		coneOS_printf("Hit MAX_ITERS, solution may be inaccurate\n"); 
	}
    coneOS_printf("Time taken: %.4f seconds\n",info->time/1e3);

	for(i = 0; i < _lineLen_; ++i) {
		coneOS_printf("-");
	}
    coneOS_printf("\n");

    if (info->stint == INFEASIBLE) {
        coneOS_printf("Certificate of primal infeasibility:\n");
        coneOS_printf("|A'y|_2 * |b|_2 = %2e\n", info->resDual);
        coneOS_printf("dist(y, K*) = 0\n");
        coneOS_printf("b'y = %.4f\n", info->dobj);
    } 
    else if (info->stint == UNBOUNDED) {
        coneOS_printf("Certificate of dual infeasibility:\n");
        coneOS_printf("|Ax + s|_2 * |c|_2 = %2e\n", info->resPri);
        coneOS_printf("dist(s, K) = 0\n");
        coneOS_printf("c'x = %.4f\n", info->pobj);
    }
    else {
        coneOS_printf("Error metrics:\n");
        coneOS_printf("|Ax + s - b|_2 / (1 + |b|_2) = %2e\n|A'y + c|_2 / (1 + |c|_2) = %2e\n",info->resPri, info->resDual);
        coneOS_printf("|c'x + b'y| / (1 + |c'x| + |b'y|) = %2e\n", info->relGap); 
        coneOS_printf("dist(s, K) = 0, dist(y, K*) = 0, s'y = 0\n");
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
