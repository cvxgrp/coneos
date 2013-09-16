void normalize(Data * d, Work * w, Cone * k){
    
	double * D = coneOS_calloc(d->m, sizeof(double));
	double * E = coneOS_calloc(d->n, sizeof(double));
	
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
	// scale the rows with D
	for(i = 0; i < d->n; ++i){
		for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
			d->Ax[j] /= D[d->Ai[j]];
		}
	}
	// calculate and scale by col norms, E
	for (i = 0; i < d->n; ++i){
		E[i] = fmax(calcNorm(&(d->Ax[d->Ap[i]]),d->Ap[i+1] - d->Ap[i]),1e-6);
		scaleArray(&(d->Ax[d->Ap[i]]), 1.0/E[i], d->Ap[i+1] - d->Ap[i]);
	}
	// scale b
	for (i = 0; i < d->m; ++i){
		d->b[i] /= D[i];
	}
	w->sc_b = 1/fmax(calcNorm(d->b,d->m),1e-6);
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
	w->sc_c = meanNormRowA/fmax(calcNorm(d->c,d->n),1e-6);
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
