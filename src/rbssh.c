#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <float.h>  	// required for DBL_EPSILON
#include <math.h>
#include "rbssh.h"

void sortN(int *ind0, double *x, int n, int dd){
	int i, j, MaxInd, d, *ind;
	double tmp;
	ind = (int*)malloc(sizeof(int)*n);
	for(i=0;i<n;i++) ind[i] = i;

	d = (dd==n?dd-1:dd);
	for(i=0;i<d;i++)
	{
		tmp = x[0]; MaxInd = ind[0];
		for(j=1;j<n-i;j++)
		{
			if(x[j]<tmp)
			{
				x[j-1] = x[j];
				x[j] = tmp;
				ind[j-1] = ind[j];
				ind[j] = MaxInd;
			}
			else
			{
				tmp = x[j];
				MaxInd = ind[j];
			}
		}	
	}
	for(j=0;j<dd;j++) ind0[j] = ind[n-j-1];
	free(ind);
}

int MatrixInvSymmetric(double *a,int n){ 
	int i,j,k,m;
    double w,g,*b;
    b = (double*)malloc(n*sizeof(double));

    for (k=0; k<=n-1; k++){ 
        w=a[0];
        if (fabs(w)+1.0==1.0){ 
            free(b); return(-2);
        }
        m=n-k-1;
        for (i=1; i<=n-1; i++){ 
            g=a[i*n]; b[i]=g/w;
            if (i<=m) b[i]=-b[i];
            for (j=1; j<=i; j++)
                a[(i-1)*n+j-1]=a[i*n+j]+g*b[j];
        }
        a[n*n-1]=1.0/w;
        for (i=1; i<=n-1; i++)
            a[(n-1)*n+i-1]=b[i];
    }
    for (i=0; i<=n-2; i++)
        for (j=i+1; j<=n-1; j++)
            a[i*n+j]=a[j*n+i];
    free(b);
    return(2);
}

int MinInd(double *x, int n){
	int i,ind=0;
	double s=x[0];
	for(i=1;i<n;i++){
		if(s>x[i]){ s=x[i]; ind=i;}
	}
	return(ind);
}

void tAbyB(double *outMatrix, const double *A, const double *B, int n, int p, int q){
    int i,j,k;
    double temp;
	for (i = 0; i<p; i++){
		for (k = 0; k<q; k++){
			temp = 0.0;
			for (j = 0; j < n; j++)
				temp += A[i + j*p] * B[k + j*q];
			outMatrix[i*q + k] = temp;
		}
	}
}

void AbyB(double *outMatrix, const double *A, const int *B, int n, int p, int q){
    int i,j,k;
    double temp;
	for (i = 0; i<n; i++){
		for (k = 0; k<q; k++){
			temp = 0.0;
			for (j = 0; j < p; j++)
				temp += A[i*p + j] * B[k + j*q];
			outMatrix[i*q + k] = temp;
		}
	}
}

int UpTriangularInv(double *B, int n, double *A){
	int i,j,k;
	for(i=0;i<n;i++)
		if(fabs(A[i*n+i])<1e-4) return(0);
	for(i=0;i<n;i++) B[i*n+i] = 1;
	for(j=1;j<n;j++)for(i=0;i<j;i++)B[j*n+i] = 0;

	for(i=n-1;i>=0;i--)//rows
	{
		if(A[i*n+i]!=1)
			for(j=i;j<n;j++)
				B[j*n+i] = B[j*n+i]/A[i*n+i];
		if(i>0)
		{
			for(j=i;j<n;j++)// columns
				for(k=0;k<i;k++)// rows
					B[j*n+k] = B[j*n+k] - A[i*n+k]*B[j*n+i];
		}
	}
	return(1);
}

void QRDecompN(double *E, double *R, double *x, int n, int p){
	double *Z, *znorm;
	double  tmp, tmp1;
	int i,j, k;
	
	Z 		= (double*)malloc(sizeof(double)*n*p);
	znorm 	= (double*)malloc(sizeof(double)*p);

	// calculate the first column
	tmp = 0;
	for(i=0;i<n;i++){
		Z[i] = x[i];
		tmp += Z[i]*Z[i];		
	}
	znorm[0] = sqrt(tmp);
	tmp = 0;
	for(i=0;i<n;i++){
		E[i] = x[i]/znorm[0];
		tmp  += E[i]*x[i];
	}
	R[0] = tmp;

	//iteration from j=1...p	
	for(j=1;j<p;j++){		
		for(k=0;k<j;k++){
			tmp=0;	for(i=0;i<n;i++) 	tmp += E[k*n+i]*x[j*n+i];
			R[j*p+k] = tmp;
		}
		tmp1 = 0;
		for(i=0;i<n;i++){
			tmp = 0; for(k=0;k<j;k++) 	tmp += R[j*p+k]*E[k*n+i];
			Z[j*n+i] = x[j*n+i] - tmp;
			tmp1 += pow(Z[j*n+i],2);
		}
		znorm[j] = sqrt(tmp1);
		tmp1 = 0;
		for(i=0;i<n;i++) 	E[j*n+i] 	= Z[j*n+i]/znorm[j];
		for(i=0;i<n;i++) 	tmp1 		+= E[j*n+i]*x[j*n+i];
		R[j*p+j] = tmp1;
	}
	free(Z); free(znorm);
}

double objectfun0(double *x, double *y, double *theta, int *delta, int n, int p, int q){
	// input:
	// x in R^{p*n}
	// y in R^{q*n}
	// theta in R^{p*q}
	// delta in R^{q}

	// output
	// obj ------ RSS
	int i,j,k;
	double obj=0.0, obj1, tmp;
	for(k=0;k<q;k++){
		if(delta[k]){
			for(i=0;i<n;i++){
				tmp = 0.0;
				for(j=0;j<p;j++)
					tmp 	+= x[j*n+i]*theta[j*q+k]; 
				obj1 	= y[k*n+i] - tmp;
				obj 	+= obj1*obj1;
			}
		}
	}
	return obj;
}

double objectfun(double *y, double *qy, int *delta, int n, int p, int q){
	// input:
	// qy in R^{p*q}
	// y in R^{q*n}
	// delta in R^{q}

	// output
	// obj ------ RSS
	int i,j;
	double obj=0.0, tmp, tmp1;
	for(j=0;j<q;j++){
		tmp = tmp1 = 0.0;
		for(i=0;i<n;i++)
			tmp		+= y[i+j*n]*y[i+j*n];
		for(i=0;i<p;i++)
			tmp1	+= qy[q*i+j]*qy[q*i+j];
		if(delta[j])
			obj 	+= tmp - tmp1;
	}
	return obj;
}

double flip1(double *Q, double *b, int *x0, int q, int maxstep){
	// input:
	// Q in R^{q*q}
	// b in R^{q}
	// x0 in R^{q}

	// output:
	// x0 in R^q
	// f0 --- the objective function
	int i,step=0,jk=0;
	double f0;
	double *Qx,*diagQ,*betak;
	Qx 		= (double*)malloc(sizeof(double)*q);    // Qx = Q*x in R^{q}
	diagQ 	= (double*)malloc(sizeof(double)*q);  	// diagQ = diag(Q) in R^{q}
	betak 	= (double*)malloc(sizeof(double)*q);  	// betak in R^{q}

	AbyB(Qx,Q,x0,q,q,1);
	f0 = 0.0;
	for(i=0;i<q;i++){
		f0 += (Qx[i]+2*b[i])*x0[i];
		diagQ[i] = Q[i*q+i];
	}

	while(step<maxstep){
		step++;
		for(i=0;i<q;i++) betak[i] = 2*(1-2*x0[i])*(Qx[i]+b[i]) + diagQ[i];
		jk = MinInd(betak,q);
		if(betak[jk]>=0.0) break;

		f0 += betak[jk];
		for(i=0;i<q;i++) Qx[i] += (1-2*x0[jk])*Q[i+jk*q];
		x0[jk] = 1 - x0[jk];
	}
	free(Qx);
	free(diagQ);
	free(betak);
	return f0;
}

double flip2(double *Q, double *b, int *x0, int q, int maxstep){
	// input:
	// Q in R^{q*q}
	// b in R^{q}
	// x0 in R^{q}

	// output:
	// x0 in R^q
	// f0 --- the objective function
	int i,j,step=0,ik=0,jk=0;
	double f0;
	double *Qx,*diagQ,*betak,*deltak;
	Qx 		= (double*)malloc(sizeof(double)*q);    // Qx = Q*x in R^{q}
	diagQ 	= (double*)malloc(sizeof(double)*q);  	// diagQ = diag(Q) in R^{q}
	betak 	= (double*)malloc(sizeof(double)*q);  	// betak in R^{q}
	deltak 	= (double*)malloc(sizeof(double)*q*q);  // betak in R^{q*q}

	AbyB(Qx,Q,x0,q,q,1);
	f0 = 0.0;
	for(i=0;i<q;i++){
		f0 += (Qx[i]+2*b[i])*x0[i];
		diagQ[i] = Q[i*q+i];
	}
	
	while(step<maxstep){
		step++;
		for(i=0;i<q;i++) betak[i] = (1-2*x0[i])*(Qx[i]+b[i]) + diagQ[i]/2;
		for(i=0;i<q;i++){
			for(j=0;j<q;j++){
				deltak[i*q+j] = betak[j] + betak[i] + (1-2*x0[i])*(1-2*x0[j])*Q[i*q+j];
			}
		}
		for(i=0;i<q;i++) deltak[i*q+i] -= diagQ[i]+betak[i];
		j = MinInd(deltak,q*q);
		ik = j%q;
		jk = (j-ik)/q;
		if(deltak[j]>=0.0) break;
		f0 += 2*deltak[j];

		if(ik==jk){
			for(i=0;i<q;i++) Qx[i] += (1-2*x0[ik])*Q[i+ik*q];
			x0[ik] = 1 - x0[ik];
		}
		else{	
			for(i=0;i<q;i++) Qx[i] += (1-2*x0[jk])*Q[i+jk*q] + (1-2*x0[ik])*Q[i+ik*q];
			x0[ik] = 1 - x0[ik];
			x0[jk] = 1 - x0[jk]; 
		}
	}
	free(Qx);
	free(diagQ);
	free(betak);
	free(deltak);
	return f0;
}

double flip12(double *Q, double *b, int *x0, int q, int maxstep){
	// input:
	// Q in R^{q*q}
	// b in R^{q}
	// x0 in R^{q}

	// output:
	// x0 in R^q
	// f0 --- the objective function
	int i,j,step=0,ik=0,jk=0;
	double f0;
	double *Qx,*diagQ,*betak,*deltak;
	Qx 		= (double*)malloc(sizeof(double)*q);    // Qx = Q*x in R^{q}
	diagQ 	= (double*)malloc(sizeof(double)*q);  	// diagQ = diag(Q) in R^{q}
	betak 	= (double*)malloc(sizeof(double)*q);  	// betak in R^{q}
	deltak 	= (double*)malloc(sizeof(double)*q*q);  // betak in R^{q*q}

	AbyB(Qx,Q,x0,q,q,1);
	f0 = 0.0;
	for(i=0;i<q;i++){
		f0 += (Qx[i]+2*b[i])*x0[i];
		diagQ[i] = Q[i*q+i];
	}
	
	while(step<maxstep){
		step++;
		for(i=0;i<q;i++) betak[i] = 2*(1-2*x0[i])*(Qx[i]+b[i]) + diagQ[i];
		jk = MinInd(betak,q);
		if(betak[jk]>=0.0) break;

		f0 += betak[jk];
		for(i=0;i<q;i++) Qx[i] += (1-2*x0[jk])*Q[i+jk*q];
		x0[jk] = 1 - x0[jk];  
	}

	step = 0;
	while(step<maxstep){
		step++;
		for(i=0;i<q;i++) betak[i] = (1-2*x0[i])*(Qx[i]+b[i]) + diagQ[i]/2;
		for(i=0;i<q;i++){
			for(j=0;j<q;j++){
				deltak[i*q+j] = betak[j] + betak[i] + (1-2*x0[i])*(1-2*x0[j])*Q[i*q+j];
			}
		}
		for(i=0;i<q;i++) deltak[i*q+i] -= diagQ[i]+betak[i];
		j = MinInd(deltak,q*q);
		ik = j%q;
		jk = (j-ik)/q;
		if(deltak[j]>=0.0) break;
		f0 += 2*deltak[j];

		if(ik==jk){
			for(i=0;i<q;i++) Qx[i] += (1-2*x0[ik])*Q[i+ik*q];
			x0[ik] = 1 - x0[ik];
		}
		else{	
			for(i=0;i<q;i++) Qx[i] += (1-2*x0[jk])*Q[i+jk*q] + (1-2*x0[ik])*Q[i+ik*q];
			x0[ik] = 1 - x0[ik];
			x0[jk] = 1 - x0[jk]; 
		}
	}
	free(Qx);
	free(diagQ);
	free(betak);
	free(deltak);
	return f0;
}

double pythag(double a, double b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

double SIGN(const double a, const double b){
	return (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

void tred2(double **z, double *d, double *e, int n, int yesvecs)
{
    int l,k,j,i;
    double scale,hh,h,g,f;
    for (i=n-1;i>0;i--) {
        l=i-1;
        h=scale=0.0;
        if (l > 0) {
            for (k=0;k<i;k++)
                scale += fabs(z[i][k]);
            if (scale == 0.0)
                e[i]=z[i][l];
            else {
                for (k=0;k<i;k++) {
                    z[i][k] /= scale;
                    h += z[i][k]*z[i][k];
                }
                f=z[i][l];
                g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i]=scale*g;
                h -= f*g;
                z[i][l]=f-g;
                f=0.0;
                for (j=0;j<i;j++) {
                    if (yesvecs)
                        z[j][i]=z[i][j]/h;
                    g=0.0;
                    for (k=0;k<j+1;k++)
                        g += z[j][k]*z[i][k];
                    for (k=j+1;k<i;k++)
                        g += z[k][j]*z[i][k];
                    e[j]=g/h;
                    f += e[j]*z[i][j];
                }
                hh=f/(h+h);
                for (j=0;j<i;j++) {
                    f=z[i][j];
                    e[j]=g=e[j]-hh*f;
                    for (k=0;k<j+1;k++)
                        z[j][k] -= (f*e[k]+g*z[i][k]);
                }
            }
        } else
            e[i]=z[i][l];
        d[i]=h;
    }
    if (yesvecs) d[0]=0.0;
    e[0]=0.0;
    for (i=0;i<n;i++) {
        if (yesvecs) {
            if (d[i] != 0.0) {
                for (j=0;j<i;j++) {
                    g=0.0;
                    for (k=0;k<i;k++)
                        g += z[i][k]*z[k][j];
                    for (k=0;k<i;k++)
                        z[k][j] -= g*z[k][i];
                }
            }
            d[i]=z[i][i];
            z[i][i]=1.0;
            for (j=0;j<i;j++) z[j][i]=z[i][j]=0.0;
        } else {
            d[i]=z[i][i];
        }
    }
}

void tqli(double **z, double *d, double *e, int n, int yesvecs)
{
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    const double EPS=DBL_EPSILON;
    for (i=1;i<n;i++) e[i-1]=e[i];
    e[n-1]=0.0;
    for (l=0;l<n;l++) {
        iter=0;
        do {
            for (m=l;m<n-1;m++) {
                dd=fabs(d[m])+fabs(d[m+1]);
                if (fabs(e[m]) <= EPS*dd) break;
            }
            if (m != l) {
                if (iter++ == 30) {
                    printf("Too many iterations in tqli");
                    exit(1);
                }
                g=(d[l+1]-d[l])/(2.0*e[l]);
                r=pythag(g,1.0);
                g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=l;i--) {
                    f=s*e[i];
                    b=c*e[i];
                    e[i+1]=(r=pythag(f,g));
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m]=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2.0*c*b;
                    d[i+1]=g+(p=s*r);
                    g=c*r-b;
                    if (yesvecs) {
                        for (k=0;k<n;k++) {
                            f=z[k][i+1];
                            z[k][i+1]=s*z[k][i]+c*f;
                            z[k][i]=c*z[k][i]-s*f;
                        }
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l]=g;
                e[m]=0.0;
            }
        } while (m != l);
    }
}

void eigencmp(double *Z, int n, double gamma)
{
	int i, j, k;
	double **z, tmp;
    z = (double**)malloc(sizeof(double*)*n);
    for (i = 0; i < n; ++i) {
        z[i] = (double*)malloc(n*sizeof(double));
        for (j = 0; j < n; ++j) z[i][j] = Z[j*n+i];
    }
    double *e = (double*)calloc(n, sizeof(double));
    double *d = (double*)calloc(n, sizeof(double));

    // svd decomposition for z; eigenvalues saved in d, left eigenvectors saved in z
    tred2(z, d, e, n, 1);
    tqli(z, d, e, n, 1);

    for (i = 0; i < n; ++i) d[i] = pow(d[i], gamma);

    for (i = 0; i < n; ++i) {
    	for (j = i; j < n; ++j) {
    		tmp = 0.0;
    		for (k = 0; k < n; ++k) tmp += z[i][k]*z[j][k]*d[k];
    		Z[j*n+i] = tmp;
    	    Z[i*n+j] = tmp;
    	}
    }
    // Z = z%*%diag(d)%*%t(z);
    for (i = 0; i < n; ++i) free(z[i]);
    free(e);
    free(d);
    free(z);
}
