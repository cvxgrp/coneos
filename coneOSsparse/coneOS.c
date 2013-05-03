/* coneos 1.0 */
#include "coneOS.h"

static int _lineLen_;
// constants and data structures
static const char* HEADER[] = {
  "Iter", 
  " norm(u - u_t) ",
  " norm(u - u_prev) ",
  " eps primal ",
  " eps dual ",
};
static const int HEADER_LEN = 5;

// to hold residual information
struct residuals {
  double resPri;
  double resDual;
  double epsPri;
  double epsDual;
};

// forward declare inline declarations
static inline void updateDualVars(Data * d, Work * w);
static inline void projectCones(Data * d,Work * w,Cone * k);
static inline void sety(Data * d, Work * w, Sol * sol);
static inline void setx(Data * d, Work * w, Sol * sol);
static inline int getSolution(Data * d, Work * w, Sol * sol, Info * info);
static inline void getInfo(Data * d, Work * w, Sol * sol, Info * info, struct residuals * r, int status);
static inline void printSummary(Data * d,Work * w,int i, struct residuals *r);
static inline void printHeader();
static inline void printFooter(Info * info); 
static inline void printSol(Data * d, Sol * sol, Info * info);
static inline void freeWork(Work * w);
static inline void projectLinSys(Data * d,Work * w);
static inline Work * initWork(Data * d, Cone * k);
static inline void unNormalize(Data *d, Work * w);
static inline void calcResiduals(Data *d, Work * w,  struct residuals * r); 


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
  struct residuals r = { -1, -1, -1, -1};
  Work * w = initWork(d,k);
  if(d->VERBOSE) {
    printHeader();
  } /* coneOS: */
  for (i=0; i < d->MAX_ITERS; ++i){
    memcpy(w->u_prev, w->u, w->l*sizeof(double));
   
   	projectLinSys(d,w);
    projectCones(d,w,k);
    updateDualVars(d,w);

    calcResiduals(d,w,&r);
    if (r.resPri < r.epsPri && r.resDual < r.epsDual) break; 
    if (d->VERBOSE && i % 100 == 0) printSummary(d,w,i, &r);
  }
  int status = getSolution(d,w,sol,info);
  info->iter = i;
  getInfo(d,w,sol,info,&r,status);
  if(d->VERBOSE) {
    printSummary(d,w,i,&r);
    printFooter(info);
    //printSol(d,sol,info);
  }
  if(d->NORMALIZE) unNormalize(d,w);
  freeWork(w);
  return status;
}

static inline void getInfo(Data * d, Work * w, Sol * sol, Info * info, struct residuals * r,int status){
	info->time = tocq();
	info->presid = r->resPri;
	info->dresid = r->resDual;
	info->pobj = innerProd(d->c,sol->x,d->n);
	info->dobj = -innerProd(d->b,sol->y,d->m);
	if (status == 2)
		info->pobj = NAN;
	else if(status == 1)
		info->dobj = NAN;
	info->gap = info->pobj-info->dobj;
}

static inline void unNormalize(Data *d, Work * w){
	int i, Anz = d->Anz;
	// unscale data
	for(i = 0; i < Anz; ++i) {
		d->Ax[i] /= (w->dual_scale)*(w->primal_scale);
	}
	for(i = 0; i < d->n; ++i) {
		d->c[i] /= w->primal_scale;
	}
	for(i = 0; i < d->m; ++i) {
		d->b[i] /= w->dual_scale;
	}
}

static inline void normalize(Data * d, Work * w){
	int i, Anz = d->Anz;
	// scale A,b,c
	double ds, ps, normA = 0.0;
	// frobenius norm
	for(i = 0; i < Anz; ++i) {
		normA = (normA > fabs(d->Ax[i])) ? normA : fabs(d->Ax[i]);
		//normA += (d->Ax[i]*d->Ax[i]);//((double)d->m*d->n);
	}
	normA = sqrt(normA);
	ds = pow((double)1.0/normA, (double)(d->n)/((double)(d->m + d->n)));
	ps = pow((double)1.0/normA, (double)(d->m)/((double)(d->m + d->n)));

	for(i = 0; i < Anz; ++i) {
		d->Ax[i] *= ds*ps;
	}
	for(i = 0; i < d->n; ++i) {
		d->c[i] *= ps;
	}
	for(i = 0; i < d->m; ++i) {
		d->b[i] *= ds;
	}
	w->dual_scale = ds;
	w->primal_scale = ps;
}

static inline Work * initWork(Data *d, Cone * k) {
  Work * w = coneOS_malloc(sizeof(Work));
  w->l = d->n+d->m+1;
  w->u = coneOS_calloc(w->l,sizeof(double));
  w->u[w->l-1] = 100;
  w->v = coneOS_calloc(w->l,sizeof(double));
  w->v[w->l-1] = 100;
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

  if(d->NORMALIZE) {
  	normalize(d,w);
  }
  else {
	  w->dual_scale = 1.0;
	  w->primal_scale = 1.0;
  }
  return w;
}

static inline void calcResiduals(Data *d ,Work *w, struct residuals * r){ 
  r->resPri = calcNormDiff(w->u,w->u_t,w->l);
  r->resDual = calcNormDiff(w->u,w->u_prev,w->l);
  r->epsPri = sqrt(w->l)*d->EPS_ABS + 
    d->EPS_REL*fmax(calcNorm(w->u,w->l),calcNorm(w->u_t,w->l));
  r->epsDual = sqrt(w->l)*d->EPS_ABS + d->EPS_REL*calcNorm(w->v,w->l);
}

static inline void projectLinSys(Data * d,Work * w){
  //memcpy(w->u_t,w->u,w->l*sizeof(double));
  //addScaledArray(w->u_t,w->v,w->l,1.0);
  //addScaledArray(w->u_t,w->h,w->l-1,-w->u_t[w->l-1]);
  //solveLinSys(d, w, w->u_t);
  //double sc = innerProd(w->u_t,w->h,w->l-1)/(w->fac+1);
  //addScaledArray(w->u_t, w->beta, w->l-1, -sc);
  //w->u_t[w->l-1] += sc;

  memcpy(w->u_t,w->u,w->l*sizeof(double));
  addScaledArray(w->u_t,w->v,w->l,1.0);
  addScaledArray(w->u_t,w->h,w->l-1,-w->u_t[w->l-1]);
  addScaledArray(w->u_t, w->h, w->l-1, -innerProd(w->u_t,w->g,w->l-1)/(w->gTh+1));
  scaleArray(&(w->u_t[d->n]),-1,d->m);
  solveLinSys(d, w, w->u_t, w->u);
  w->u_t[w->l-1] += innerProd(w->u_t,w->h,w->l-1);
}

static inline void freeWork(Work * w){
  freePriv(w);
  if(w->Xs) coneOS_free(w->Xs);
  if(w->Z) coneOS_free(w->Z);
  if(w->e) coneOS_free(w->e);
  if(w->u) coneOS_free(w->u);
  if(w->v) coneOS_free(w->v);
  if(w->u_t) coneOS_free(w->u_t);
  if(w->u_prev) coneOS_free(w->u_prev);
  if(w->h) coneOS_free(w->h);
  if(w->g) coneOS_free(w->g);
  if(w) coneOS_free(w);
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
  for(i = 0; i < w->l; ++i) { 
    w->v[i] += (w->u[i] - d->ALPH*w->u_t[i] - (1.0 - d->ALPH)*w->u_prev[i]); 
  }
}

void projectCones(Data *d,Work * w,Cone * k){
  int i;
  for(i = 0; i < w->l; ++i){
    w->u[i] += d->ALPH*(w->u_t[i] - w->u[i]) - w->v[i];
  }
  /* u = [x;y;tau] */
  projCone(&(w->u[d->n]),k,w);
  if (w->u[w->l-1]<0.0) w->u[w->l-1] = 0.0;
}

static inline int getSolution(Data * d,Work * w,Sol * sol, Info * info){
  double tau = (w->u[w->l-1]+w->u_t[w->l-1])/2;
  double kap = w->v[w->l-1];
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
  //memcpy(sol->y, w->z + w->yi, d->m*sizeof(double));
  int i;
  for(i = 0; i < d->m; ++i) {
    sol->y[i] = 0.5 * w->dual_scale * (w->u[i + d->n]+w->u_t[i + d->n]);
  }
}

static inline void setx(Data * d,Work * w, Sol * sol){
  sol->x = coneOS_malloc(sizeof(double)*d->n);
  //memcpy(sol->x, w->z, d->n*sizeof(double));
  int i;
  for(i = 0; i < d->n; ++i) {
    sol->x[i] = w->primal_scale * w->u[i];
  }
}

static inline void printSummary(Data * d,Work * w,int i, struct residuals *r){
  // coneOS_printf("Iteration %i, primal residual %4f, primal tolerance %4f\n",i,err,EPS_PRI);
  coneOS_printf("%*i | ", (int)strlen(HEADER[0]), i);
  coneOS_printf("%*.4f   ", (int)strlen(HEADER[1]), r->resPri); // p_res
  coneOS_printf("%*.4f   ", (int)strlen(HEADER[2]), r->resDual); // d_res
  coneOS_printf("%*.4f   ", (int)strlen(HEADER[3]), r->epsPri); // full(p_inf));
  coneOS_printf("%*.4f\n", (int)strlen(HEADER[4]), r->epsDual);//full(d_inf));
}

static inline void printHeader() {
  int i;  
  _lineLen_ = 0;
  for(i = 0; i < HEADER_LEN - 1; ++i) {
    _lineLen_ += strlen(HEADER[i]) + 3;
  }
  _lineLen_ += strlen(HEADER[HEADER_LEN-1]);
  
  for(i = 0; i < _lineLen_; ++i) {
    coneOS_printf("-");
  }
  coneOS_printf("\nconeOS 1.0\n");
  for(i = 0; i < _lineLen_; ++i) {
    coneOS_printf("-");
  }
  coneOS_printf("\n");
  for(i = 0; i < HEADER_LEN - 1; ++i) {
    coneOS_printf("%s | ", HEADER[i]);
  }
  coneOS_printf("%s\n", HEADER[HEADER_LEN-1]);
  for(i = 0; i < _lineLen_; ++i) {
    coneOS_printf("=");
  }
  coneOS_printf("\n");
}

static inline void printFooter(Info * info) {
 int i;  
  for(i = 0; i < _lineLen_; ++i) {
    coneOS_printf("-");
  }
  coneOS_printf("\nStatus: %s\n",info->status); 
  coneOS_printf("Primal objective value: %4f\n",info->pobj);
  coneOS_printf("Dual objective value: %4f\n",info->dobj);
  coneOS_printf("Duality gap: %e\n", info->gap);
  coneOS_printf("Time taken %4f ms\n",info->time);

  for(i = 0; i < _lineLen_; ++i) {
    coneOS_printf("=");
  }
  coneOS_printf("\n");
}
