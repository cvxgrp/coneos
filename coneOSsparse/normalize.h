#ifndef NORMAL_H_GUARD
#define NORMAL_H_GUARD

#define MIN_SCALE 1e-2
#define MAX_SCALE 1e3

void normalize(Data * d, Work * w, Cone * k){
    
	double * D = coneOS_calloc(d->m, sizeof(double));
	double * E = coneOS_calloc(d->n, sizeof(double));
	
	int i, j, count;
	double wrk;
	
	// heuristic rescaling, seems to do well with a scaling of about 4
	w->scale =  4; //fmax( fmin( sqrt( d->n * d->m / d->Anz ) , MAX_SCALE) , 1);

	// calculate row norms
	for(i = 0; i < d->n; ++i){
		for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
			wrk = d->Ax[j];
			D[d->Ai[j]] += wrk*wrk;
		}
	}
	for (i=0; i < d->m; ++i){
		D[i] = sqrt(D[i]); // just the norms
	}
    // mean of norms of rows across each cone 
    count = k->l+k->f;
	for(i = 0; i < k->qsize; ++i)
    {
		wrk = 0;
	    /*
        for (j = count; j < count + k->q[i]; ++j){
        	wrk = fmax(wrk,D[j]);
		}
        */
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
        /*
        for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
        	wrk = fmax(wrk,D[j]);
		}
        */
		for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
        	wrk += D[j]*D[j];
		}
        wrk = sqrt(wrk);
        wrk /= k->s[i];
        for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
        	D[j] = wrk;
		}
		count += (k->s[i])*(k->s[i]);
    }
    
    for(i = 0; i < k->ep + k->ed; ++i)
    {
        wrk = D[count]/3 + D[count + 1]/3 + D[count + 2]/3;
        D[count] = wrk;
        D[count + 1] = wrk;
        D[count + 2] = wrk;
        count += 3;
    }

    for (i=0; i<d->m; ++i){
        if (D[i] < MIN_SCALE) D[i] = MIN_SCALE;
        else if (D[i] > MAX_SCALE) D[i] = MAX_SCALE;

    }
	// scale the rows with D
	for(i = 0; i < d->n; ++i){
		for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
			d->Ax[j] /= D[d->Ai[j]];
		}
	}
	// calculate and scale by col norms, E
	for (i = 0; i < d->n; ++i){
		E[i] = calcNorm(&(d->Ax[d->Ap[i]]),d->Ap[i+1] - d->Ap[i]);
        if (E[i] < MIN_SCALE) E[i] = MIN_SCALE;
		else if (E[i] > MAX_SCALE) E[i] = MAX_SCALE;
        scaleArray(&(d->Ax[d->Ap[i]]), 1.0/E[i], d->Ap[i+1] - d->Ap[i]);
	}
	// scale b
	for (i = 0; i < d->m; ++i){
		d->b[i] /= D[i];
	}
	w->sc_b = 1/fmax(calcNorm(d->b,d->m),MIN_SCALE);
	scaleArray(d->b, w->sc_b, d->m);
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
	w->sc_c = meanNormRowA/fmax(calcNorm(d->c,d->n),MIN_SCALE);
	scaleArray(d->c, w->sc_c, d->n);

	w->D = D;
	w->E = E;

	coneOS_free(nms);

    // heuristic scaling factor
    scaleArray(d->Ax,w->scale,d->Anz);
    scaleArray(d->b,w->scale,d->m);
    scaleArray(d->c,w->scale,d->n);
    /*
       coneOS_printf("norm D is %4f\n", calcNorm(D,d->m));
       coneOS_printf("norm E is %4f\n", calcNorm(E,d->n));
       coneOS_printf("norm A is %4f\n", calcNorm(d->Ax,d->Anz));
       coneOS_printf("norm b is %4f\n", calcNorm(d->b,d->m));
       coneOS_printf("norm c is %4f\n", calcNorm(d->c,d->n));
     */
}

void calcScaledResids(Data * d, Work * w, struct residuals * r) {
    double * D = w->D;
    double * E = w->E;
    double * u = w->u;
    double * u_t = w->u_t;
    double * u_prev = w->u_prev;
    double tmp;
    int i;

    r->resPri = 0;
    for (i = 0; i < d->n; ++i){
        tmp = (u[i] - u_t[i])/(E[i] * w->sc_b);
        r->resPri += tmp * tmp;        
    }
    for (i = 0; i < d->m; ++i){
        tmp = (u[i + d->n] - u_t[i + d->n])/(D[i] * w->sc_c);
        r->resPri += tmp * tmp;
    }
    tmp = u[w->l-1] - u_t[w->l-1];
    r->resPri += tmp * tmp;
    r->resPri = sqrt(r->resPri);

    r->resDual = 0;
    for (i = 0; i < d->n; ++i){
        tmp = (u[i] - u_prev[i]) * E[i] / w->sc_b;
        r->resDual += tmp * tmp;        
    }
    for (i = 0; i < d->m; ++i){
        tmp = (u[i + d->n] - u_prev[i + d->n]) * D[i] / w->sc_c;
        r->resDual += tmp * tmp;
    }
    tmp = u[w->l-1] - u_t[w->l-1];
    r->resDual += tmp * tmp;
    r->resDual = sqrt(r->resDual);
}

double calcScaledNormInf(const double *a, const double * s, int l){ 
    double tmp, max = 0.0;
    int i;
    for ( i=0; i<l; ++i){
        tmp = fabs(s[i] * a[i]);
        if(tmp > max) max = tmp;
    }   
    return max;
}

void unNormalize(Data *d, Work * w, Sol * sol){
	int i, j;
	double * D = w->D;
	double * E = w->E;
	for (i = 0; i < d->n; ++i){
		sol->x[i] /= (E[i] * w->sc_b);
	}
	for (i = 0; i < d->m; ++i){
		sol->y[i] /= (D[i] * w->sc_c);
	}
	for (i = 0; i < d->m; ++i){
		sol->s[i] *= D[i]/(w->sc_b * w->scale);
	}
    for (i = 0; i < d->n; ++i){
		d->c[i] *= E[i]/(w->sc_c * w->scale);
	}
	for (i = 0; i < d->m; ++i){
		d->b[i] *= D[i]/(w->sc_b * w->scale);
	}
    for (i = 0; i < d->n; ++i){
        scaleArray(&(d->Ax[d->Ap[i]]), E[i], d->Ap[i+1] - d->Ap[i]);
    }   
    for(i = 0; i < d->n; ++i){
        for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
            d->Ax[j] *= D[d->Ai[j]];
        }   
    }
    scaleArray(d->Ax,1.0/w->scale,d->Anz);
}

#endif
