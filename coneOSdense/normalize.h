#ifndef NORMAL_H_GUARD
#define NORMAL_H_GUARD

void normalize(Data * d, Work * w, Cone * k){
    
	double * D = coneOS_calloc(d->m, sizeof(double));
	double * E = coneOS_calloc(d->n, sizeof(double));
	
	int i, j, count;
	double wrk;
	
	// heuristic rescaling, seems to do well with a scaling of about 4
	w->scale = 4.0;

    // calculate row norms
    for(i = 0; i < d->m; ++i){
        D[i] = cblas_dnrm2(d->n,&(d->Ax[i]),d->m);
    }
    // mean of norms of rows across each cone	
    count = k->l+k->f;
    for(i = 0; i < k->qsize; ++i)
    {
        wrk = 0;
        /*
        for (j = count; j < count + k->q[i]; ++j){
            wrk = fmax(wrk, D[j]);
        }   
        */
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
        /*
        for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
            wrk = fmax(wrk, D[j]);
        }
        */
        for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
            wrk += D[j]*D[j];
        }
        wrk = sqrt(wrk);
        wrk /= (k->s[i]);
        for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
            D[j] = wrk;
        }
        count += (k->s[i])*(k->s[i]);
    }
    for (i=0; i<d->m; ++i){
        if (D[i] < 1e-1) D[i] = 1;
    }
	// scale the rows with D
	for(i = 0; i < d->m; ++i){
		cblas_dscal(d->n,1.0/D[i],&(d->Ax[i]),d->m);
	}
	// calculate and scale by col norms, E
	for (i = 0; i < d->n; ++i){
		E[i] = cblas_dnrm2(d->m,&(d->Ax[i*d->m]),1);
        if (E[i] < 1e-3) E[i] = 1;
		cblas_dscal(d->m,1.0/E[i],&(d->Ax[i*d->m]),1);	
	}
	// scale b
	for (i = 0; i < d->m; ++i){
		d->b[i] /= D[i];
	}
	w->sc_b = 1/fmax(calcNorm(d->b,d->m),1e-1);
	scaleArray(d->b, w->sc_b, d->m);
	// scale c
	for (i = 0; i < d->n; ++i){
		d->c[i] /= E[i];
	}

	double meanNormRowA = 0.0;
	for(i = 0; i < d->m; ++i){
		meanNormRowA += cblas_dnrm2(d->n,&(d->Ax[i]),d->m)/d->m;
	}
	w->sc_c = meanNormRowA/fmax(calcNorm(d->c,d->n),1e-1);
	scaleArray(d->c, w->sc_c, d->n);

	w->D = D;
	w->E = E;

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

inline double calcScaledNormDiff(double * a, double * b, Data * d, Work * w) {
    double nmDiff = 0.0, tmp;
    double * D = w->D;
    double * E = w->E;
    int i;
    for (i = 0; i < d->n; ++i){
        tmp = (a[i] - b[i])/(E[i] * w->sc_b);
        nmDiff += tmp * tmp;    
    }    
    for (i = 0; i < d->m; ++i){
        tmp = (a[i + d->n] - b[i + d->n])/(D[i] * w->sc_c);
        nmDiff += tmp * tmp;
    }   
    tmp = a[w->l-1] - b[w->l-1];
    nmDiff += tmp * tmp;
    return sqrt(nmDiff);
}

void unNormalize(Data *d, Work * w, Sol * sol){
	int i;
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
 		cblas_dscal(d->m,E[i],&(d->Ax[i*d->m]),1);	
	}   
	for(i = 0; i < d->m; ++i){
		cblas_dscal(d->n,D[i],&(d->Ax[i]),d->m);
	}
	scaleArray(d->Ax,1.0/w->scale,d->Anz);
}

#endif
