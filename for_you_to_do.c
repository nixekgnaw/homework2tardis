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
    /* add your code here */
    int i, j, t, k, l, maxind, max;
    double temps;
    for (i = 0; i < n - 1; i++) {
        maxind = i;
        max = abs(A[i * n + i]);
        for (t = i + 1; t < n; t++) {
            if (abs(A[t * n + i]) > max) {
                maxind = t;
                max = abs(A[t * n + i]);
            }
        }
        if (max == 0) {
            return -1;
        } else {
            if (maxind != i) {
                temps = ipiv[i];
                ipiv[i] = ipiv[maxind];
                ipiv[maxind] = temps;
                for (l = 0; l < n; l++) {
                    temps = A[i * n + l];
                    A[i * n + l] = A[maxind * n + l];
                    A[maxind * n + l] = temps;
                }
            }
            for (j = i + 1; j < n; j++) {
                A[j * n + i] = A[j * n + i] / A[i * n + i];
                for (k = i + 1; k < n; k++) {
                    A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
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
    if (UPLO == 'L') {
        for (i = 1; i < n; i++) {
            for (j = 0; j < i; j++)
                B[ipiv[i]] -= B[ipiv[j]] * A[i * n + j];
        }
    } else if (UPLO == 'U') {
        B[ipiv[n - 1]] /= A[(n - 1) * n + n - 1];
        for (i = n - 2; i >= 0; i--)
            for (j = i + 1; j < n; j++)
                B[ipiv[i]] -= B[ipiv[j]] * A[i * n + j];
        B[ipiv[i]] /= A[i * n + i];
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
    return;
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
    return 0;
}

