#include <math.h> 		// required for sqrt(), fabs();
#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <R.h>
#include <Rinternals.h> // required for SEXP et.al.;
#include "rbssh.h"

double Flip_single(int *delta, double *theta, double *x, double *y, double *V, double ga,
						double lambda, int n, int p, int q, int isV, int isflip1){
	// input:
	// x in R^{p*n}
	// y in R^{q*n}
	// V in R^{q*q}

	// output:
	// theta in R^{p*q}
	// delta in R^{q}
	int i,j,k;
	double tmp,tmp1, obj;
	double *Q, *R, *qy, *invR,*A,*b;
	Q 		= (double*)malloc(sizeof(double)*n*p);  // Q in R^{p*n}
	R 		= (double*)malloc(sizeof(double)*p*p);  // R in R^{p*p}
	qy 		= (double*)malloc(sizeof(double)*p*q);  // qy in R^{p*q}
	invR 	= (double*)malloc(sizeof(double)*p*p);	// invR in R^{p*p}
	A 		= (double*)malloc(sizeof(double)*q*q);  // A in R^{q*q}
	b 		= (double*)malloc(sizeof(double)*q);    // b in R^{q}

	QRDecompN(Q, R, x, n, p);

	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			tmp = 0.0;
			for(i=0;i<n;i++) 	tmp += Q[j*n+i]*y[k*n+i];
			qy[q*j+k] = tmp;
		}
	}
	if(isV){
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				tmp = tmp1 = 0.0;
				for(i=0;i<n;i++)	tmp += y[k*n+i]*y[j*n+i];
				for(i=0;i<p;i++)	tmp1 += qy[i*q+k]*qy[i*q+j];
				A[k*q+j] = V[k*q+j]*(tmp - tmp1);
				if(j==k) b[k] 		= -lambda*V[k*q+k]*pow(tmp1,ga);
			}
		}
	}
	else{
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				tmp = tmp1 = 0.0;
				for(i=0;i<n;i++)	tmp += y[k*n+i]*y[j*n+i];
				for(i=0;i<p;i++)	tmp1 += qy[i*q+k]*qy[i*q+j];
				V[k*q+j] = (tmp - tmp1)/(n-p);
				A[k*q+j] = tmp - tmp1;
				if(j==k) b[k] = -lambda*pow(tmp1,ga);
			}
		}
		MatrixInvSymmetric(V,q);
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				A[k*q+j] *= V[k*q+j];
				if(j==k) b[k] *= V[k*q+k];
			}
		}
	}

	for(j=0;j<q;j++) {
		delta[j] 	= 0;
		A[j*q+j] 	-= b[j];
	}
	if(isflip1==1)			obj = flip1(A,b,delta,q,q);
	else if(isflip1==2)		obj = flip2(A,b,delta,q,q);
	else 					obj = flip12(A,b,delta,q,q);

	for(i=0;i<p*q;i++) 	theta[i] = 0.0;
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
	free(A);
	free(b);
	return obj;
}

int Flip_bic(int *delta, double *theta, double *bic, double *x, double *y, double *V, double ga,
				double *lambda, int n, int p, int q, int nlam, double tau, int isV, int isflip1, int criteria)
{
	// input:
	// x in R^{p*n}
	// y in R^{q*n}
	// v in R^{q*q}
	// lambda in R^{nlam}

	// output:
	// bic in R^{nlam}
	// theta in R^{p*q}
	// delta in R^{q*nlam}

	int i,j,k,minid,df,*deltak;
	double tmp,tmp1,minbic,obj;
	double *Q, *R, *invR, *qy, *A0, *A, *b, *bk;
	Q 		= (double*)	malloc(sizeof(double)*n*p);  	// Q in R^{p*n}
	R 		= (double*)	malloc(sizeof(double)*p*p);  	// R in R^{p*p}
	invR 	= (double*)	malloc(sizeof(double)*p*p);		// R in R^{p*p}
	qy 		= (double*)	malloc(sizeof(double)*p*q); 	// qy in R^{p*q}
	deltak 	= (int*)	malloc(sizeof(int)*q); 			// qy in R^{p*nlam}
	A0 		= (double*)	malloc(sizeof(double)*q*q);  	// A in R^{q*q}
	A 		= (double*)	malloc(sizeof(double)*q*q);  	// A in R^{q*q}
	b 		= (double*)	malloc(sizeof(double)*q);    	// b in R^{q}
	bk 		= (double*)	malloc(sizeof(double)*q);    	// bk in R^{q}

	QRDecompN(Q, R, x, n, p);

	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			tmp = 0.0;
			for(i=0;i<n;i++) 	tmp += Q[j*n+i]*y[k*n+i];
			qy[q*j+k] = tmp;
		}
	}
	if(isV){
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				tmp = tmp1 = 0.0;
				for(i=0;i<n;i++)	tmp += y[k*n+i]*y[j*n+i];
				for(i=0;i<p;i++)	tmp1 += qy[i*q+k]*qy[i*q+j];
				A0[k*q+j] 	= V[k*q+j]*(tmp - tmp1);
				if(j==k)    b[k] = V[k*q+k]*pow(tmp1,ga);
			}
		}
	}
	else{
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				tmp = tmp1 = 0.0;
				for(i=0;i<n;i++)	tmp += y[k*n+i]*y[j*n+i];
				for(i=0;i<p;i++)	tmp1 += qy[i*q+k]*qy[i*q+j];
				V[k*q+j] = (tmp - tmp1)/(n-p);
				A0[k*q+j] = tmp - tmp1;
				if(j==k) 	b[k] = pow(tmp1,ga);
			}
		}
		MatrixInvSymmetric(V,q);
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				A0[k*q+j] *= V[k*q+j];
				if(j==k) b[k] *= V[k*q+k];
			}
		}
	}

	for(k=0;k<nlam;k++){
		df = 0;
		for(i=0;i<q;i++){
			deltak[i] 	= 0;
			bk[i] 		= -lambda[k]*b[i];
			A[i*q+i] 	= A0[i*q+i] - bk[i];
		}
		if(isflip1==1)			obj = flip1(A,bk,deltak,q,q);
		else if(isflip1==2)		obj = flip2(A,bk,deltak,q,q);
		else 					obj = flip12(A,bk,deltak,q,q);
		for(j=0;j<q;j++){
			delta[k*q+j] = deltak[j];
			if(deltak[j]) 	df++;
		}
		obj = objectfun(y, qy, deltak, n, p, q);

		switch (criteria){
		case 1: bic[k] = log(obj/n/q) + 2*pow(1.0*df,tau)/n/q; 	break; 	// AIC
		case 2: bic[k] = log(obj/n/q) + log(1.0*n*q)*df/n/q; 	break;	// BIC
		case 3: bic[k] = obj*(n*q)/(n*q-df)/(n*q-df); 			break;	// GCV
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
	for(j=0;j<q;j++)	deltak[j] = delta[minid*q+j];

	for(i=0;i<p*q;i++) 	theta[i] = 0.0;
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
	free(A0);
	free(A);
	free(b);
	free(bk);
	return minid;
}

void Flip_cv(double *bic, double *x, double *y, double *xt, double *yt, double *V, double ga,
				double *lambda, int n, int nt, int p, int q, int nlam, int isV, int isflip1)
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
	double tmp,tmp1;
	double *Q, *R, *invR, *qy, *A0, *A, *b, *bk, *theta;
	Q 		= (double*)	malloc(sizeof(double)*n*p);  	// Q in R^{p*n}
	R 		= (double*)	malloc(sizeof(double)*p*p);  	// R in R^{p*p}
	invR 	= (double*)	malloc(sizeof(double)*p*p);		// R in R^{p*p}
	qy 		= (double*)	malloc(sizeof(double)*p*q); 	// qy in R^{p*q}
	deltak 	= (int*)	malloc(sizeof(int)*q); 			// qy in R^{p*nlam}
	A 		= (double*)	malloc(sizeof(double)*q*q);  	// A in R^{q*q}
	A0 		= (double*)	malloc(sizeof(double)*q*q);  	// A in R^{q*q}
	b 		= (double*)	malloc(sizeof(double)*q);    	// b in R^{q}
	bk 		= (double*)	malloc(sizeof(double)*q);    	// bk in R^{q}
	theta	= (double*)	malloc(sizeof(double)*p*q);  	// R in R^{p*p}

	QRDecompN(Q, R, x, n, p);
	for(i=0;i<p*q;i++) 	theta[i] = 0.0;
	UpTriangularInv(invR, p, R);

	for(k=0;k<q;k++){
		for(j=0;j<p;j++){
			tmp = 0.0;
			for(i=0;i<n;i++) 	tmp 	+= Q[j*n+i]*y[k*n+i];
			qy[q*j+k] = tmp;
		}
	}
	if(isV){
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				tmp = tmp1 = 0.0;
				for(i=0;i<n;i++)	tmp 	+= y[k*n+i]*y[j*n+i];
				for(i=0;i<p;i++)	tmp1 	+= qy[i*q+k]*qy[i*q+j];
				A0[k*q+j] = V[k*q+j]*(tmp - tmp1);
				if(j==k) 	b[k] = V[k*q+k]*pow(tmp1,ga);
			}
		}
	}
	else{
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				tmp = tmp1 = 0.0;
				for(i=0;i<n;i++)	tmp 	+= y[k*n+i]*y[j*n+i];
				for(i=0;i<p;i++)	tmp1 	+= qy[i*q+k]*qy[i*q+j];
				V[k*q+j] = (tmp - tmp1)/(n-p);
				A0[k*q+j] = tmp - tmp1;
				if(j==k) 	b[k] = pow(tmp1,ga);
			}
		}
		MatrixInvSymmetric(V,q);
		for(k=0;k<q;k++){
			for(j=0;j<q;j++){
				A0[k*q+j] *= V[k*q+j];
				if(j==k) b[k] *= V[k*q+k];
			}
		}
	}

	for(s=0;s<nlam;s++){
		for(i=0;i<q;i++){
			deltak[i] 	= 0;
			bk[i] 		= -lambda[s]*b[i];
			A[i*q+i] 	= A0[i*q+i] - bk[i];
		}
		if(isflip1==1)			flip1(A, bk,deltak,q,q);
		else if(isflip1==2)		flip2(A, bk,deltak,q,q);
		else 					flip12(A,bk,deltak,q,q);

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
	free(invR);
	free(qy);
	free(deltak);
	free(A0);
	free(A);
	free(b);
	free(bk);
	free(theta);
}

SEXP RBSS_FLIP(SEXP X_, SEXP Y_, SEXP V_, SEXP DIM_, SEXP PARAM_)
{
	// dimensions
	int *dim      	= INTEGER(DIM_);
	int n         	= dim[0];
	int p     	  	= dim[1];
	int q         	= dim[2];
	int isV			= dim[3];
	int isflip1   	= dim[4];
	double *param 	= REAL(PARAM_);
	double ga     	= param[0];
	double lambda 	= param[1];
	int i;

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);
	double *V		= REAL(V_);

	// Outcome
	SEXP rDELTA, rTHETA, rRSS, list, list_names;
  	PROTECT(rDELTA 	= allocVector(INTSXP, q));
  	PROTECT(rTHETA 	= allocVector(REALSXP, p*q));
  	PROTECT(rRSS   	= allocVector(REALSXP, 1));

	REAL(rRSS)[0] = Flip_single(INTEGER(rDELTA), REAL(rTHETA), x, y, V, ga, lambda, n, p, q, isV, isflip1);

	char *names[3] = {"delta", "theta", "rss"};
	PROTECT(list_names = allocVector(STRSXP, 3));
	for(i = 0; i < 3; i++)
		SET_STRING_ELT(list_names, i,  mkChar(names[i]));
	PROTECT(list = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(list, 0, rDELTA);
	SET_VECTOR_ELT(list, 1, rTHETA);
	SET_VECTOR_ELT(list, 2, rRSS);
	setAttrib(list, R_NamesSymbol, list_names);

	UNPROTECT(5);
	return list;
}

SEXP RBSS_FLIP_BIC(SEXP X_, SEXP Y_, SEXP V_, SEXP LAMBDA_, SEXP DIM_, SEXP PARAM_)
{
	// dimensions
	int *dims  		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1];
	int q     		= dims[2];
	int nlam  		= dims[3];
	int isV			= dims[4];
	int isflip1   	= dims[5];
	int criteria	= dims[6];

	double *param 	= REAL(PARAM_);
	double ga  		= param[0];
	double tau 		= param[1];
	double *lambda 	= REAL(LAMBDA_);
	int i;

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);
	double *V		= REAL(V_);

	// Outcome
	SEXP rDELTA, rTHETA, rBIC, rDF, list, list_names;
  	PROTECT(rDELTA 	= allocVector(INTSXP, q*nlam));
  	PROTECT(rTHETA 	= allocVector(REALSXP, p*q));
  	PROTECT(rBIC   	= allocVector(REALSXP, nlam));
	PROTECT(rDF   	= allocVector(INTSXP, 1));

	INTEGER(rDF)[0] = Flip_bic(INTEGER(rDELTA), REAL(rTHETA), REAL(rBIC), x, y, V, ga, lambda,
					n, p, q, nlam, tau, isV, isflip1, criteria);

	char *names[4] = {"delta", "theta", "bic", "selected"};
	PROTECT(list_names = allocVector(STRSXP, 4));
	for(i = 0; i < 4; i++)
		SET_STRING_ELT(list_names, i,  mkChar(names[i]));
	PROTECT(list = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(list, 0, rDELTA);
	SET_VECTOR_ELT(list, 1, rTHETA);
	SET_VECTOR_ELT(list, 2, rBIC);
	SET_VECTOR_ELT(list, 3, rDF);
	setAttrib(list, R_NamesSymbol, list_names);

	UNPROTECT(6);
	return list;
}

SEXP RBSS_FLIP_CV(SEXP X_, SEXP Y_, SEXP Xt_, SEXP Yt_, SEXP V_, SEXP LAMBDA_, SEXP DIM_, SEXP PARAM_)
{
	// dimensions
	int *dims  		= INTEGER(DIM_);
	int n     		= dims[0];
	int p     		= dims[1];
	int q     		= dims[2];
	int nlam  		= dims[3];
	int isV			= dims[4];
	int isflip1   	= dims[5];
	int nt			= dims[6];

	double *param 	= REAL(PARAM_);
	double ga  		= param[0];
	double *lambda 	= REAL(LAMBDA_);

	// Pointers
	double *x 		= REAL(X_);
	double *y  		= REAL(Y_);
	double *xt 		= REAL(Xt_);
	double *yt  	= REAL(Yt_);
	double *V		= REAL(V_);

	// Outcome
	SEXP rBIC, list, list_names;
  	PROTECT(rBIC   		= allocVector(REALSXP, 	nlam));
	PROTECT(list 		= allocVector(VECSXP, 	1));
	PROTECT(list_names 	= allocVector(STRSXP, 	1));

	Flip_cv(REAL(rBIC), x, y, xt, yt, V, ga, lambda, n, nt, p, q, nlam, isV, isflip1);

	SET_STRING_ELT(list_names, 	0,  mkChar("rss"));
	SET_VECTOR_ELT(list, 		0, 	rBIC);
	setAttrib(list, R_NamesSymbol, 	list_names);

	UNPROTECT(3);
	return list;
}

SEXP FLIP(SEXP A_, SEXP b_, SEXP x0_, SEXP PARAM_){
	int *param  	= INTEGER(PARAM_);
	int q			= param[0];
	int isflip1		= param[1];


	double *A 		= REAL(A_);
	double *b  		= REAL(b_);
	int *x0			= INTEGER(x0_);

	int i;
	double 	obj;

	if(isflip1==1)			obj = flip1(A,b,x0,q,q);
	else if(isflip1==2)		obj = flip2(A,b,x0,q,q);
	else 					obj = flip12(A,b,x0,q,q);

	// Outcome
	SEXP rXhat, rOBJ, list, list_names;
  	PROTECT(rXhat 	= allocVector(INTSXP, q));
  	PROTECT(rOBJ 	= allocVector(REALSXP, 1));

	for(i = 0; i < q; i++) 	INTEGER(rXhat)[i] = x0[i];
	REAL(rOBJ)[0] = obj;

	char *names[2] = {"xhat", "obj"};
	PROTECT(list_names = allocVector(STRSXP, 2));
	for(i = 0; i < 2; i++)	SET_STRING_ELT(list_names, i,  mkChar(names[i]));
	PROTECT(list = allocVector(VECSXP, 2)); 
	SET_VECTOR_ELT(list, 0, rXhat);
	SET_VECTOR_ELT(list, 1, rOBJ);
	setAttrib(list, R_NamesSymbol, list_names); 

	UNPROTECT(4);
	return list;
}

