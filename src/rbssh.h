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

double objectfun(double *y, double *qy, int *delta, int n, int p, int q);

double flip1(double *Q, double *b, int *x0, int q, int maxstep);

double flip2(double *Q, double *b, int *x0, int q, int maxstep);

double flip12(double *Q, double *b, int *x0, int q, int maxstep);

#endif // RBSS_H_INCLUDED
