#include <math.h> 		// required for sqrt(), fabs();
#include <stdio.h>  	// required for exit
#include <stdlib.h> 	// required for malloc(), free();
#include <string.h> 	// required for memcpy()
#include <R.h>
#include <Rinternals.h> // required for SEXP et.al.;
#include "rbssh.h"

// Check if the element is in the array
int is_in_array(int *arr, int size, int element) {
    for (int i = 0; i < size; i++) {
        if (arr[i] == element) {
            return 1;
        }
    }
    return 0;
}

// omega_k
double compute_omega_k(double *Y, double *X, int n, int p, int k) {
    int i, j;
    double omega_k = 0.0, sum = 0.0;
    for (j = 0; j < n; j++) {
        sum = 0.0;
        for (i = 0; i < n; i++) {
            if (Y[i] < Y[j]) {
                sum += X[i * p + k];
            }
        }
        omega_k += pow(sum / n, 2);
    }
    return (n) / (double)((n - 1) * (n - 2)) * omega_k;
}

// SIRS
void SIRS(int *ind, double *Y, double *X, int n, int p, int N, int d, int dd) {
    int i, k, l, soft_size, ind_size;
    double C_d;
    double *omega   = (double*)malloc(sizeof(double) * p);
    int *ind_hard   = (int*)malloc(sizeof(int) * N);
    int *ind_soft   = (int*)malloc(sizeof(int) * p);
    int *ind_whole  = (int*)malloc(sizeof(int) * p);
    int *ind1       = (int*)malloc(sizeof(int) * 2);
    double *z       = (double*)malloc(sizeof(double) * n * d);
    double *omega_z = (double*)malloc(sizeof(double) * d);

    // Initialize the ind array to -1 (indicating unassigned)
    for (i = 0; i < p; i++) {
        ind[i] = -1;
    }

    for (k = 0; k < p; k++) {
        omega[k] = compute_omega_k(Y, X, n, p, k);
    }

    // hard threshold
    sortN(ind_hard, omega, p, N); // sort

    // soft threshold
    GetRNGstate();  // Obtain R random number seeds
    for (l = 0; l < d; l++) {
        for (i = 0; i < n; i++) {
            z[i * d + l] = unif_rand();  // Generate random numbers
        }
        omega_z[l] = compute_omega_k(Y, z, n, d, l);
    }
    PutRNGstate();  // Release R random number seeds

    sortN(ind1, omega_z, d, 2); // sort
    C_d = omega_z[ind1[0]];

    sortN(ind_whole, omega, p, p);
    soft_size = 0;
    for (i = 0; i < p; i++) {
        k = ind_whole[i];
        if (omega[k] > C_d) {
            ind_soft[soft_size++] = k;
        }
    }

    ind_size = 0;
    for (i = 0; i < N; i++) {
        if (!is_in_array(ind, ind_size, ind_hard[i])) {
            ind[ind_size++] = ind_hard[i];
        }
    }
    for (i = 0; i < soft_size; i++) {
        if (!is_in_array(ind, ind_size, ind_soft[i])) {
            ind[ind_size++] = ind_soft[i];
        }
    }

    free(omega);
    free(ind_hard);
    free(ind_soft);
    free(ind_whole);
    free(ind1);
    free(z);
    free(omega_z);
}

SEXP _SIRS(SEXP Y_, SEXP X_, SEXP DIM_, SEXP N_, SEXP d_, SEXP dd_) {
    int *dims = INTEGER(DIM_);
    int n = dims[0];
    int p = dims[1];
    int N = INTEGER(N_)[0];
    int d = INTEGER(d_)[0];
    int dd = INTEGER(dd_)[0]; // The final output feature count
    int i;

    double *y = REAL(Y_);
    double *x = REAL(X_);

    double *omega = (double *)malloc(p * sizeof(double));
    for (i = 0; i < p; i++) {
        omega[i] = compute_omega_k(y, x, n, p, i);
    }

    SEXP rA;
    PROTECT(rA = allocVector(INTSXP, p));
    int *A = INTEGER(rA); // Get pointer
    int A_size = 0; // The actual number of selected features

    SIRS(A, y, x, n, p, N, d, dd);

    for (i = 0; i < p; i++) {
        if (A[i] != -1) {
            A_size++;
        }
    }

    // If the actual number of selected features is less than dd, only output the selected features
    SEXP rA_final;
    if (A_size < dd) {
        PROTECT(rA_final = allocVector(INTSXP, A_size)); // Assign an integer vector of length A_size
        memcpy(INTEGER(rA_final), A, A_size * sizeof(int)); // Copy the actual selected features
    } else {
        PROTECT(rA_final = allocVector(INTSXP, dd));
        memcpy(INTEGER(rA_final), A, dd * sizeof(int));
    }

    SEXP rOmega;
    PROTECT(rOmega = allocVector(REALSXP, p));
    memcpy(REAL(rOmega), omega, p * sizeof(double));

    SEXP list, list_names;
    char *names[2] = {"omega", "indn"};
    PROTECT(list_names = allocVector(STRSXP, 2));
    for (i = 0; i < 2; i++) {
        SET_STRING_ELT(list_names, i, mkChar(names[i]));
    }
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, rOmega);
    SET_VECTOR_ELT(list, 1, rA_final);
    setAttrib(list, R_NamesSymbol, list_names);

    UNPROTECT(5);
    return list;
}