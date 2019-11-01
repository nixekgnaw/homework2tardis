#include "../include/for_you_to_do.h"

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf(double *A, int *ipiv, int n) {
    int i, j, t, k, l, maxind, tempi;
    double temps, max;
    for (i = 0; i < n - 1; i++) {
        maxind = i;
        max = fabs(A[i * n + i]);
        for (t = i + 1; t < n; t++) {
            if (fabs(A[t * n + i]) > max) {
                maxind = t;
                max = fabs(A[t * n + i]);
            }
        }
        if (max == 0) {
            return -1;
        } else {
            if (maxind != i) {
                tempi = ipiv[i];
                ipiv[i] = ipiv[maxind];
                ipiv[maxind] = tempi;
                for (l = 0; l < n; l++) {
                    temps = A[i * n + l];
                    A[i * n + l] = A[maxind * n + l];
                    A[maxind * n + l] = temps;
                }
            }
            for (j = i + 1; j < n; j++) {
                A[j * n + i] /= A[i * n + i];
                for (k = i + 1; k < n; k++) {
                    A[j * n + k] -= A[j * n + i] * A[i * n + k];
                }
            }
        }
    }
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv) {

    /* add your code here */
    int i, j;
    double sum;
    if (UPLO == 'L') {
        double *y;
        y = (double *) malloc(sizeof(double) * n);
        y[0] = B[ipiv[0]];
        for (i = 1; i < n; i++) {
            sum = 0;
            for (j = 0; j < i; j++)
                sum += y[j] * A[i * n + j];
            y[i] = B[ipiv[i]] - sum;
        }
        for (i = 0; i < n; i++)
            B[i] = y[i];
    } else if (UPLO == 'U') {
        B[n - 1] /= A[(n - 1) * n + n - 1];
        for (i = n - 2; i >= 0; i--) {
            sum = 0;
            for (j = i + 1; j < n; j++)
                sum += B[j] * A[i * n + j];
            B[i] = (B[i] - sum) / A[i * n + i];
        }
    }
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b) {
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    int ii, jj, kk, i1, j1, k1;
    for (ii = i; ii < i + b && ii < n; ii += b)
        for (jj = j; jj < j + b && jj < n; jj += b)
            for (kk = k; kk < k + b && kk < n; kk += b)
                for (i1 = ii; i1 < ii + b && i1 < n; i1 += 3) {

                    for (j1 = jj; j1 < jj + b && j1 < n; j1 += 3) {

                        register int t = i1 * n + j1;
                        register int tt = t + n;
                        register int ttt = tt + n;
                        register double c00 = C[t];
                        register double c01 = C[t + 1];
                        register double c02 = C[t + 2];
                        register double c10 = C[tt];
                        register double c11 = C[tt + 1];
                        register double c12 = C[tt + 2];
                        register double c20 = C[ttt];
                        register double c21 = C[ttt + 1];
                        register double c22 = C[ttt + 2];
                        for (k1 = kk; k1 < kk + b && k1 < n; k1 += 3) {
                            register int ta = i1 * n + k1;
                            register int tta = ta + n;
                            register int ttta = tta + n;
                            register int tb = k1 * n + j1;
                            register int ttb = tb + n;
                            register int tttb = ttb + n;

                            register double r1 = A[ta];     // a00
                            register double r2 = A[tta];    // a10
                            register double r3 = A[ttta];   // a20
                            register double r4 = B[tb];     // b00
                            register double r5 = B[tb + 1]; // b01
                            register double r6 = B[tb + 2]; // b02

                            c00 -= r1 * r4;
                            c01 -= r1 * r5;
                            c02 -= r1 * r6;
                            c10 -= r2 * r4;
                            c11 -= r2 * r5;
                            c12 -= r2 * r6;
                            c20 -= r3 * r4;
                            c21 -= r3 * r5;
                            c22 -= r3 * r6;
                            r1 = A[ta + 1];
                            r2 = A[tta + 1];
                            r3 = A[ttta + 1];
                            r4 = B[ttb];
                            r5 = B[ttb + 1];
                            r6 = B[ttb + 2];
                            c00 -= r1 * r4;
                            c01 -= r1 * r5;
                            c02 -= r1 * r6;
                            c10 -= r2 * r4;
                            c11 -= r2 * r5;
                            c12 -= r2 * r6;
                            c20 -= r3 * r4;
                            c21 -= r3 * r5;
                            c22 -= r3 * r6;
                            r1 = A[ta + 2];
                            r2 = A[tta + 2];
                            r3 = A[ttta + 2];
                            r4 = B[tttb];
                            r5 = B[tttb + 1];
                            r6 = B[tttb + 2];
                            c00 -= r1 * r4;
                            c01 -= r1 * r5;
                            c02 -= r1 * r6;
                            c10 -= r2 * r4;
                            c11 -= r2 * r5;
                            c12 -= r2 * r6;
                            c20 -= r3 * r4;
                            c21 -= r3 * r5;
                            c22 -= r3 * r6;
                        }
                        C[t] = c00;
                        C[t + 1] = c01;
                        C[t + 2] = c02;
                        C[tt] = c10;
                        C[tt + 1] = c11;
                        C[tt + 2] = c12;
                        C[ttt] = c20;
                        C[ttt + 1] = c21;
                        C[ttt + 2] = c22;
                    }
                }
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) {
    int ib, i, t, k, j, l, maxind, tempi;
    double max, sum, temps;
    for (ib = 0; ib < n; ib += b) {
        for (i = ib; i < ib + b && i < n; i++) {
            max = fabs(A[i * n + i]);
            maxind = i;
            for (t = i + 1; t < n; t++) {
                if (fabs(A[t * n + i]) > max) {
                    maxind = t;
                    max = fabs(A[t * n + i]);
                }
            }

            if (max == 0) {
                return -1;
            } else {
                if (maxind != i) {
                    tempi = ipiv[i];
                    ipiv[i] = ipiv[maxind];
                    ipiv[maxind] = tempi;

                    for (j = 0; j < n; j++) {
                        temps = A[i * n + j];
                        A[i * n + j] = A[maxind * n + j];
                        A[maxind * n + j] = temps;
                    }
                }
            }
            for (k = i + 1; k < n; k++) {
                A[k * n + i] = A[k * n + i] / A[i * n + i];
                for (l = i + 1; l < ib + b && l < n; l++)
                    A[k * n + l] -= A[i * n + l] * A[k * n + i];
            }
        }

        for (i = ib; i < ib + b && i < n; i++) {
            for (j = ib + b; j < n; j++) {
                sum = 0;
                for (k = ib; k < i; k++) {
                    sum += A[i * n + k] * A[k * n + j];
                }
                A[i * n + j] = (A[i * n + j]-sum)/A[i*n+i];
            }
        }

        mydgemm(A, A, A, n, i, j, ib, b);
    }

    return 0;
}

