#ifndef RBSS_H_INCLUDED
#define RBSS_H_INCLUDED

#define FMAX(a,b) ((a) > (b) ? (a) : (b))

void sortN(int *ind0, double *x, int n, int dd);

int MatrixInvSymmetric(double *a,int n);

int MinInd(double *x, int n);

void tAbyB(double *outMatrix, const double *A, const double *B, int n, int p, int q);

void AbyB(double *outMatrix, const double *A, const int *B, int n, int p, int q);

int UpTriangularInv(double *B, int n, double *A);

void QRDecompN(double *E, double *R, double *x, int n, int p);

double objectfun0(double *x, double *y, double *theta, int *delta, int n, int p, int q);
	// input:
	// x in R^{p*n}
	// y in R^{q*n}
	// theta in R^{p*q}
	// delta in R^{q}

double objectfun(double *y, double *qy, int *delta, int n, int p, int q);
	// input:
	// qy in R^{p*q}
	// y in R^{q*n}
	// delta in R^{q}


#endif // RBSS_H_INCLUDED
