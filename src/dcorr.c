#include <math.h> 		// required for sqrt(), fabs();
#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <R.h>
#include <Rinternals.h> // required for SEXP et.al.;
#include "rbssh.h"

void DistCorrGroup(double *dcov, double *y, double *x, int py, int px, int n){
	int i,j,k;
	double xx1, yy1, sxy1_, sxy2_, sxy3_, sxx=0, syy=0, S1, S2, S22, S3;
	double *sxy1, *sxy2, *sxy3;

	sxy1 = (double*)malloc(sizeof(double)*n);
	sxy2 = (double*)malloc(sizeof(double)*n);
	sxy3 = (double*)malloc(sizeof(double)*n);

	for(i=0;i<n;i++){
		sxy1_ = sxy2_ = sxy3_ = 0.0;
		for(j=0;j<n;j++){
			yy1 = xx1 = 0; 
			for(k=0;k<px;k++) xx1 += (x[k*n+j] - x[k*n+i])*(x[k*n+j] - x[k*n+i]);
			for(k=0;k<py;k++) yy1 += (y[k*n+j] - y[k*n+i])*(y[k*n+j] - y[k*n+i]);
			sxx += xx1;
			syy += yy1;

			yy1 = sqrt(yy1); 
			xx1 = sqrt(xx1);
			sxy1_ += xx1*yy1;
			sxy2_ += xx1;
			sxy3_ += yy1;
		}
		sxy1[i] = sxy1_/n;
		sxy2[i] = sxy2_/n;
		sxy3[i] = sxy3_/n;
	}
	sxx = sxx/n/n; 
	syy = syy/n/n;

	// calculate the dcovXY//
	S1 = S2 = S3 = S22 = 0.0;
	for(i=0;i<n;i++){
		S1 += sxy1[i];
		S2 += sxy2[i]; 	S22 += sxy3[i];
		S3 += sxy2[i]*sxy3[i];
	}
	dcov[3] = sqrt((S1 + S22*S2/n - 2*S3)/n);

	// calculate the dvarXX//
	S1 = sxx; 
	S3=0;
	for(i=0;i<n;i++)	S3 += sxy2[i]*sxy2[i];
	dcov[1] = sqrt(S1 + S2*S2/n/n - 2*S3/n);

	// calculate the dvarYY//
	S1 = syy; S3=0;
	for(i=0;i<n;i++)	S3 += sxy3[i]*sxy3[i];
	dcov[2] = sqrt(S1 + S22*S22/n/n - 2*S3/n);
	dcov[0] = dcov[3]/sqrt(dcov[1]*dcov[2]);
	free(sxy1);
	free(sxy2);
	free(sxy3);
}

void DistCorr(int *ind0, double *corr, double *y, double *x, int py, int p, int n, int dd){
	int i,j;
	double *dcov, *xj;
	dcov = (double*)malloc(sizeof(double)*4);
	xj = (double*)malloc(sizeof(double)*n);

	for(j=0;j<p;j++){
		for(i=0;i<n;i++) xj[i] = x[j*n+i];
		DistCorrGroup(dcov, y, xj, py, 1, n);
		corr[j] = dcov[0];
	}
	sortN(ind0, corr, p, dd);

	free(xj);
	free(dcov);
}

SEXP DCORR_SCREEN(SEXP X_, SEXP Y_, SEXP DIM_, SEXP d_)
{
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1]; 
	int q     		= dims[2];
	int d			= INTEGER(d_)[0];
	int i;

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);

	// Outcome
	SEXP rDCOV, rNTOP, list, list_names;
	PROTECT(rDCOV 	= allocVector(REALSXP, p));
	PROTECT(rNTOP 	= allocVector(INTSXP, d));
	DistCorr(INTEGER(rNTOP), REAL(rDCOV), y, x, q, p, n, d);

	char *names[4] = {"dcor", "indn"};
	PROTECT(list_names = allocVector(STRSXP, 2));
	for(i = 0; i < 2; i++)
		SET_STRING_ELT(list_names, i,  mkChar(names[i]));
	PROTECT(list = allocVector(VECSXP, 2)); 
	SET_VECTOR_ELT(list, 0, rDCOV);
	SET_VECTOR_ELT(list, 1, rNTOP);
	setAttrib(list, R_NamesSymbol, list_names); 

	UNPROTECT(4);
	return list;
}

SEXP DCORR(SEXP X_, SEXP Y_, SEXP DIM_)
{
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1]; 
	int q     		= dims[2];

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);

	// Outcome
	SEXP rDCOV;
  	PROTECT(rDCOV 	= allocVector(REALSXP, 4));
	DistCorrGroup(REAL(rDCOV), y, x, q, p, n);
	UNPROTECT(1);
	return rDCOV;
}



