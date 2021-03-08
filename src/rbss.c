#include <math.h> 		// required for sqrt(), fabs();
#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <R.h>
#include <Rinternals.h> // required for SEXP et.al.;
#include "rbssh.h"

void PvaluesMarginal(double *pvalues, double *sigma2, double *x, double *y, int n, int p, int q){
	// input:
	// x in R^{p*n}
	// y in R^{q*n}

	// output:
	// pvalues in R^{q}
	// sigma2 in R^{q}
	
	int i,j,k;
	double tmp,tmp1,s;
	double *Q, *R, *qy;
	Q 		= (double*)malloc(sizeof(double)*n*p);  // Q in R^{p*n}
	R 		= (double*)malloc(sizeof(double)*p*p);  // R in R^{p*p}
	qy 		= (double*)malloc(sizeof(double)*p*q);  // qy in R^{p*q}

	QRDecompN(Q, R, x, n, p);
	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			tmp = 0.0; 
			for(i=0;i<n;i++) 	tmp += Q[j*n+i]*y[k*n+i];
			qy[q*j+k] = tmp;
		}
	}

	for(j=0;j<q;j++){
		tmp = tmp1 = 0.0;
		for(i=0;i<n;i++)
			tmp 	+= y[i+j*n]*y[i+j*n];
		for(i=0;i<p;i++)
			tmp1 	+= qy[q*i+j]*qy[q*i+j];
		s 	= (tmp-tmp1)/(n-p);
		pvalues[j] 	= tmp1/s;
		sigma2[j] 	= s;
	}	
	free(Q);
	free(R);
	free(qy);
}

double Delta(int *delta, double *qy, double *y, double ga, double lambda, int n, int q, int p){
	// input:
	// y in R^{q*n}
	// qy in R^{p*q}

	// output:
	// delta in R^q
	// obj --- the objective function
	int i,j;
	double tmp,tmp1,obj=0.0;
	for(j=0;j<q;j++){
		tmp = tmp1 = 0.0;
		for(i=0;i<n;i++)	tmp 	+= y[i+j*n]*y[i+j*n];
		for(i=0;i<p;i++)	tmp1 	+= qy[q*i+j]*qy[q*i+j];
		if(tmp-tmp1<=lambda*pow(tmp1,ga)){
			delta[j] 	= 1;
			obj 		+= tmp - tmp1;
		}
		else{	
			delta[j] 	= 0;
			obj 		+= tmp;
		}
	}
	return obj;
}

double Delta_single(int *delta, double *theta, double *x, double *y, double ga, double lambda, int n, int p, int q){
	// input:
	// x in R^{p*n}
	// y in R^{q*n}

	// output:
	// theta in R^{p*q}
	// delta in R^{q}	
	int i,j,k;
	double tmp, obj;
	double *Q, *R, *qy, *invR;// *iden, *Xh;
	Q 		= (double*)malloc(sizeof(double)*n*p);   // Q in R^{p*n}
	R 		= (double*)malloc(sizeof(double)*p*p);   // R in R^{p*p}
	qy 		= (double*)malloc(sizeof(double)*p*q);   // qy in R^{p*q}
	invR 	= (double*)malloc(sizeof(double)*p*p);	 // R in R^{p*p}

	QRDecompN(Q, R, x, n, p); // x = R*Q, where R is LowTriangular
	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			tmp = 0.0; 
			for(i=0;i<n;i++) 	tmp += Q[j*n+i]*y[k*n+i];
			qy[q*j+k] = tmp;
		}
	}
	obj = Delta(delta, qy, y, ga, lambda, n, q, p);

	for(i=0;i<p*q;i++) theta[i] = 0.0;
	UpTriangularInv(invR, p, R);
	for(k=0;k<q;k++){
		if(delta[k]){
			for(j=0;j<p;j++){
				tmp = 0.0;
				for(i=j;i<p;i++)	tmp += invR[i*p+j]*qy[i*q+k];
				theta[j*q+k] = tmp;
			}
		}
	}	
	free(Q);
	free(R);
	free(invR);
	free(qy);
	return obj;
}

int Delta_bic(int *delta, double *theta, double *bic, double *x, double *y, double ga, 
				double *lambda, int n, int p, int q, int nlam, double tau, int criteria)
{
	// input:
	// x in R^{p*n}
	// y in R^{q*n}
	// lambda in R^{nlam}

	// output:
	// bic in R^{nlam}
	// theta in R^{p*q}
	// delta in R^{nlam*q}	

	int i,j,k,minid,df,*deltak;
	double tmp,minbic;
	double *Q, *R, *invR, *qy;
	Q 		= (double*)	malloc(sizeof(double)*n*p);  	// Q in R^{p*n}
	R 		= (double*)	malloc(sizeof(double)*p*p);  	// R in R^{p*p}
	invR 	= (double*)	malloc(sizeof(double)*p*p);  	// R in R^{p*p}
	qy 		= (double*)	malloc(sizeof(double)*p*q);  	// qy in R^{p*q}
	deltak 	= (int*)	malloc(sizeof(int)*q);     		// qy in R^{p}
	
	QRDecompN(Q, R, x, n, p);	// x = R*Q, where R is LowTriangular

	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			tmp = 0.0; 
			for(i=0;i<n;i++) 	tmp += Q[j*n+i]*y[k*n+i];
			qy[q*j+k] = tmp;
		}
	}

	
	for(k=0;k<nlam;k++){
		df = 0;
		tmp = Delta(deltak, qy, y, ga, lambda[k], n, q, p);
		for(j=0;j<q;j++){ 
			delta[k*q+j] = deltak[j];
			if(deltak[j]) df++;
		}
		switch (criteria){
		case 1: bic[k] = log(tmp/n/q) + 2*pow(1.0*(p+1)*df,tau)/n/q; 	break; 	// AIC
		case 2: bic[k] = log(tmp/n/q) + log(n*q)*df*(p+1)/n/q; 			break;	// BIC
		case 3: bic[k] = tmp*(n*q)/(n*q-df)/(n*q-df); 					break;	// GCV
		}
		
	}
	minbic = bic[0];
	minid = 0;
	for(k=1;k<nlam;k++){
		if(bic[k]<minbic){ 
			minbic = bic[k];
			minid = k;
		}
	}
	for(j=0;j<q;j++)	deltak[j] 	= delta[minid*q+j];

	for(i=0;i<p*q;i++) 	theta[i] 	= 0.0;
	UpTriangularInv(invR, p, R);
	for(k=0;k<q;k++){
		if(deltak[k]){
			for(j=0;j<p;j++){
				tmp = 0.0;
				for(i=j;i<p;i++)	tmp += invR[i*p+j]*qy[i*q+k];
				theta[j*q+k] = tmp;
			}
		}
	}
	free(Q);
	free(R);
	free(invR);
	free(qy);
	free(deltak);
	return minid;
}

void Delta_cv(double *bic, double *x, double *y, double *xt, double *yt, double ga, double *lambda, int n, int nt, int p, int q, int nlam)
{
	// input:
	// x in R^{p*n} ------- training x 
	// y in R^{q*n} ------- training y
	// xt in R^{p*n} ------ test x
	// yt in R^{q*n} ------ test y
	// lambda in R^{nlam}

	// output:
	// bic in R^{nlam}	

	int i,j,k,s,*deltak;
	double tmp;
	double *Q, *R, *qy, *invR, *theta;
	Q 		= (double*)	malloc(sizeof(double)*n*p);  	// Q in R^{p*n}
	R 		= (double*)	malloc(sizeof(double)*p*p);  	// R in R^{p*p}
	qy 		= (double*)	malloc(sizeof(double)*p*q);  	// qy in R^{p*q}
	invR 	= (double*)	malloc(sizeof(double)*p*p);  	// R in R^{p*p}
	deltak 	= (int*)	malloc(sizeof(int)*q);     	 	// qy in R^{p}
	theta	= (double*)	malloc(sizeof(double)*p*q);  	// R in R^{p*p}
	
	QRDecompN(Q, R, x, n, p);	// x = R*Q, where R is LowTriangular
	for(i=0;i<p*q;i++) 	theta[i] = 0.0;
	UpTriangularInv(invR, p, R);

	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			tmp = 0.0; 
			for(i=0;i<n;i++) 	tmp += Q[j*n+i]*y[k*n+i];
			qy[q*j+k] = tmp;
		}
	}
	
	for(s=0;s<nlam;s++){
		Delta(deltak, qy, y, ga, lambda[s], n, q, p);	
		for(k=0;k<q;k++){
			if(deltak[k]){
				for(j=0;j<p;j++){
					tmp = 0.0;
					for(i=j;i<p;i++)	tmp += invR[i*p+j]*qy[i*q+k];
					theta[j*q+k] = tmp;
				}
			}
		}
		bic[s] = objectfun0(xt, yt, theta, deltak, nt, p, q);
	}

	free(Q);
	free(R);
	free(qy);
	free(invR);
	free(deltak);
	free(theta);
}

SEXP RBSS(SEXP X_, SEXP Y_, SEXP DIM_, SEXP PARAM_)
{
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1]; 
	int q     		= dims[2];
	double *param 	= REAL(PARAM_);
	double ga 		= param[0];
	double lambda 	= param[1];
	int i;

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);

	// Outcome
	SEXP rDELTA, rTHETA, rRSS, list, list_names;
  	PROTECT(rDELTA 		= allocVector(INTSXP, 	q));
  	PROTECT(rTHETA 		= allocVector(REALSXP, 	p*q));
  	PROTECT(rRSS   		= allocVector(REALSXP, 	1));
	PROTECT(list 		= allocVector(VECSXP, 	3));
	PROTECT(list_names 	= allocVector(STRSXP, 	3)); 

	REAL(rRSS)[0] = Delta_single(INTEGER(rDELTA), REAL(rTHETA), x, y, ga, lambda, n, p, q);

	char *names[3] = {"delta", "theta", "rss"};	
	for(i = 0; i < 3; i++)
		SET_STRING_ELT(list_names, i,  mkChar(names[i]));
	
	SET_VECTOR_ELT(list, 0, rDELTA);
	SET_VECTOR_ELT(list, 1, rTHETA);
	SET_VECTOR_ELT(list, 2, rRSS);  
	setAttrib(list, R_NamesSymbol, list_names); 

	UNPROTECT(5);
	return list;
}

SEXP RBSS_BIC(SEXP X_, SEXP Y_, SEXP LAMBDA_, SEXP DIM_, SEXP PARAM_)
{
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1]; 
	int q     		= dims[2];
	int nlam  		= dims[3];
	int criteria	= dims[4];

	double *param 	= REAL(PARAM_);
	double ga  		= param[0];
	double tau 		= param[1];
	double *lambda 	= REAL(LAMBDA_);
	int i;

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);

	// Outcome
	SEXP rDELTA, rTHETA, rBIC, rDF, list, list_names;
  	PROTECT(rDELTA 		= allocVector(INTSXP, q*nlam));
  	PROTECT(rTHETA 		= allocVector(REALSXP, p*q));
  	PROTECT(rBIC   		= allocVector(REALSXP, nlam));
	PROTECT(rDF   		= allocVector(INTSXP, 1));
	PROTECT(list 		= allocVector(VECSXP, 4)); 
	PROTECT(list_names 	= allocVector(STRSXP, 4));

	INTEGER(rDF)[0] = Delta_bic(INTEGER(rDELTA), REAL(rTHETA), REAL(rBIC), x, y, ga, lambda, n, p, q, nlam, tau, criteria);

	

	char *names[4] 		= {"delta", "theta", "bic", "selected"};
	for(i = 0; i < 4; i++)
		SET_STRING_ELT(list_names, i,  mkChar(names[i]));
	
	SET_VECTOR_ELT(list, 0, rDELTA);
	SET_VECTOR_ELT(list, 1, rTHETA);
	SET_VECTOR_ELT(list, 2, rBIC); 
	SET_VECTOR_ELT(list, 3, rDF); 
	setAttrib(list, R_NamesSymbol, list_names); 

	UNPROTECT(6);
	return list;
}

SEXP RBSS_CV(SEXP X_, SEXP Y_, SEXP Xt_, SEXP Yt_, SEXP LAMBDA_, SEXP DIM_, SEXP PARAM_)
{
	// dimensions
	int *dims 		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1]; 
	int q     		= dims[2];
	int nlam  		= dims[3];
	int nt			= dims[4];

	double *param 	= REAL(PARAM_);
	double ga  		= param[0];
	double *lambda 	= REAL(LAMBDA_);

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);
	double *xt 		= REAL(Xt_);
	double *yt  	= REAL(Yt_);
	
	// Outcome
	SEXP rBIC, list, list_names;
  	PROTECT(rBIC   		= allocVector(REALSXP, 	nlam));
	PROTECT(list 		= allocVector(VECSXP, 	1)); 
	PROTECT(list_names 	= allocVector(STRSXP, 	1));

	Delta_cv(REAL(rBIC), x, y, xt, yt, ga, lambda, n, nt, p, q, nlam);
	
	SET_STRING_ELT(list_names, 	0,  mkChar("rss"));	
	SET_VECTOR_ELT(list, 		0, 	rBIC); 
	setAttrib(list, R_NamesSymbol, 	list_names); 

	UNPROTECT(3);
	return list;
}

SEXP PVALUES(SEXP X_, SEXP Y_, SEXP DIM_){
	int *dim      	= INTEGER(DIM_);
	int n         	= dim[0];
	int p     	  	= dim[1]; 
	int q         	= dim[2];

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);

	// Outcome
	SEXP rPvalues, rSigma2, list, list_names;
  	PROTECT(rPvalues 	= allocVector(REALSXP, 	q));
	PROTECT(rSigma2 	= allocVector(REALSXP, 	q));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));
	PROTECT(list 		= allocVector(VECSXP, 	2)); 

	PvaluesMarginal(REAL(rPvalues), REAL(rSigma2), x, y, n, p, q);

	SET_STRING_ELT(list_names, 	0,	mkChar("Tn"));
	SET_STRING_ELT(list_names, 	1,  mkChar("Sigma2"));
	SET_VECTOR_ELT(list, 		0, 	rPvalues);
	SET_VECTOR_ELT(list, 		1, 	rSigma2);
	setAttrib(list, R_NamesSymbol, 	list_names); 

	UNPROTECT(4);
	return list;
}



