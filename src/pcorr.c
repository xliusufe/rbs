#include <math.h> 		// required for sqrt(), fabs();
#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <R.h>
#include <Rinternals.h> // required for SEXP et.al.;
#include "rbssh.h"

void PearsonCorrGroup(double* corr, double* y, double* x, int p, int n) {
    int i, j;
    double* mean_x = (double*)calloc(p, sizeof(double));
    double* var_x = (double*)calloc(p, sizeof(double));
    double* cov_xy = (double*)calloc(p, sizeof(double));
    double mean_y, var_y;
    mean_y = var_y = 0.0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) mean_x[j] += x[i * p + j];
        mean_y += y[i];
    }
    for (j = 0; j < p; j++) mean_x[j] /= n;
    mean_y /= n;

    for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
            var_x[j] += (x[i * p + j] - mean_x[j]) * (x[i * p + j] - mean_x[j]);
        }
        var_y += (y[i] - mean_y) * (y[i] - mean_y);
        for (j = 0; j < p; j++) {
            cov_xy[j] += (x[i * p + j] - mean_x[j]) * (y[i] - mean_y);
        }
    }

    // Pearson correlation
    for (j = 0; j < p; j++) {
        if (var_x[j] > 0 && var_y > 0) {
            corr[j] = cov_xy[j] / (sqrt(var_x[j]) * sqrt(var_y));
        }
        else {
            corr[j] = 0.0; // If the variance is 0 and the correlation is 0
        }
    }

    free(mean_x);
    free(var_x);
    free(cov_xy);
}

void PearsonCorrRank(int *ind0, double* corr, double* y, double* x, int p, int n, int dd){
    int j;
    double *Fabs_Corr;
    Fabs_Corr = (double*)malloc(sizeof(double)*p);

    PearsonCorrGroup(corr, y, x, p, n);
    for (j = 0; j < p; j++) {
        Fabs_Corr[j] = fabs(corr[j]);
    }
    sortN(ind0, Fabs_Corr, p, dd); // sort

    free(Fabs_Corr);
}

SEXP PEARSON_SCREEN(SEXP Y_, SEXP X_, SEXP DIM_, SEXP d_) {
    int *dims = INTEGER(DIM_);
    int n = dims[0];
    int p = dims[1];
    int d = INTEGER(d_)[0];
    int i;

    double *y = REAL(Y_);
    double *x = REAL(X_);

    SEXP rCORR, rNTOP, list, list_names;
    PROTECT(rCORR = allocVector(REALSXP, p));
    PROTECT(rNTOP = allocVector(INTSXP, d));
    PearsonCorrRank(INTEGER(rNTOP), REAL(rCORR), y, x, p, n, d);

    char *names[2] = {"pcor", "indn"};
    PROTECT(list_names = allocVector(STRSXP, 2));
    for (i = 0; i < 2; i++) {
        SET_STRING_ELT(list_names, i, mkChar(names[i]));
    }
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, rCORR);
    SET_VECTOR_ELT(list, 1, rNTOP);
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(4);
    return list;
}

SEXP PEARSON_CORR(SEXP Y_, SEXP X_, SEXP DIM_)
{
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1];

	// Pointers
	double *y       = REAL(Y_);
    double *x       = REAL(X_);

	// Outcome
	SEXP rCORR;
  	PROTECT(rCORR = allocVector(REALSXP, p));
	PearsonCorrGroup(REAL(rCORR), y, x, p, n);
	UNPROTECT(1);
	return rCORR;
}