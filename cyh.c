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
int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    // printf("%s/n","I'm in mydgetrf");
    double max;
    int i, maxind, t, j, k;
    for (i=0; i<n; i++) {
        // ensure A can use LU decomposition
        maxind = i;
        max = fabs(A[i*n+i]);
        // printf("初始 max 的值为 %e", max);
        for (t=i+1; t<n; t++)
            if (fabs(A[t*n+i])>max) {
                maxind = t;
                max = fabs(A[t*n+i]);
                // printf("max 的值为 %e", max);
            }
        if (max == 0) {
            printf("%s/n","LU factoration failed: coefficient matirx is singular..？");
            return -1;
        }
        else {
            if (maxind != i) {
                int temps, index;
                double *tempv=(double *) malloc(sizeof(double) * n);
                temps = ipiv[i];
                ipiv[i] = ipiv[maxind];
                ipiv[maxind] = temps;
                // swap rows
                for (index = 0;index < n; index++){
                    tempv[index] = A[i*n+index];
                    A[i*n+index] = A[maxind*n+index];
                    A[maxind*n+index] = tempv[index];
                } 
            }
        }
        // factorization
        for (j = i+1; j<n; j++) {
            A[j*n+i] = A[j*n+i]/A[i*n+i];
            for (k=i+1; k<n; k++)
                A[j*n+k] -= A[j*n+i]*A[i*n+k];
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
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    // forward substitution
    int i, j;
    double sum;
    double *y=(double *) malloc(sizeof(double) * n);
    if (UPLO == 'L' ) {
        y[0] = B[ipiv[0]];
        for (i=1; i<n; i++) {
            sum = 0.0;
            for (j= 0; j<i; j++)
                sum += y[j]*A[i*n+j];
            y[i] = B[ipiv[i]] - sum;
        }
        int k;
        for(k=0;k<n;k++)
            B[k]=y[k];
    }
    // backward subsititution
    
    if (UPLO == 'U') {
        double *x=(double *) malloc(sizeof(double) * n);
        x[n-1] = B[n-1] / A[(n-1)*n+n-1];
        for (i=n-2; i>=0; i--) {
            sum = 0.0;
            for (j=i+1; j<n; j++)
                sum += x[j]*A[i*n+j];
            x[i] = (B[i]-sum)/A[i*n+i];
        }
        int k;
        for(k=0;k<n;k++)
            B[k]=x[k];
    }
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    int ii,kk,jj;
    int i1,j1,k1;
    for (ii=i; ii<i+b && ii<n; ii+=b)
        for (jj=j; jj<j+b && jj<n; jj+=b)
            for (kk=k; kk<k+b && kk<n; kk+=b)
                for (i1=ii; i1<ii+b && i1<n; i1+=3)
                    for (j1=jj; j1<jj+b && j1<n; j1+=3) {
                        register int t= i1*n+j1; // [i,j]
                        register int tt= t+n; // [i+1,j]
                        register int ttt= tt+n; // [i+2,j]
                        register double c00 = C[t]; register double c01 = C[t+1]; register double c02 = C[t+2];
                        register double c10 = C[tt]; register double c11 = C[tt+1]; register double c12 = C[tt+2];
                        register double c20 = C[ttt]; register double c21 = C[ttt+1]; register double c22 = C[ttt+2];

                        for (k1=kk; k1<k+b && k1<n; k1+=3) {
                            /* 3 by 3 mini matrix multiplication using registers */
                            register int ta = i1*n+k1; // [i,k]
                            register int tta = ta+n; // [i+1,k]
                            register int ttta = tta+n; // [i+2,k]
                            register int tb = k1*n+j1; // [k,j]
                            register int ttb = tb+n; // [k+1,j]
                            register int tttb = ttb+n; // [k+2,j]

                            register double a00 = A[ta]; register double a10 = A[tta]; register double a20 = A[ttta];
                            register double b00 = B[tb]; register double b01 = B[tb+1]; register double b02 = B[tb+2];

                            c00 -= a00*b00; c01 -= a00*b01; c02 -= a00*b02;
                            c10 -= a10*b00; c11 -= a10*b01; c12 -= a10*b02;
                            c20 -= a20*b00; c21 -= a20*b01; c22 -= a20*b02;

                            a00 = A[ta+1]; a10 = A[tta+1]; a20 = A[ttta+1];
                            b00 = B[ttb]; b01 = B[ttb+1]; b02 = B[ttb+2];

                            c00 -= a00*b00; c01 -= a00*b01; c02 -= a00*b02;
                            c10 -= a10*b00; c11 -= a10*b01; c12 -= a10*b02;
                            c20 -= a20*b00; c21 -= a20*b01; c22 -= a20*b02;

                            a00 = A[ta+2]; a10 = A[tta+2]; a20 = A[ttta+2];
                            b00 = B[tttb]; b01 = B[tttb+1]; b02 = B[tttb+2];

                            c00 -= a00*b00; c01 -= a00*b01; c02 -= a00*b02;
                            c10 -= a10*b00; c11 -= a10*b01; c12 -= a10*b02;
                            c20 -= a20*b00; c21 -= a20*b01; c22 -= a20*b02;
                        }
                        C[t] = c00; C[t+1] = c01; C[t+2] = c02;
                        C[tt] = c10; C[tt+1] = c11; C[tt+2] = c12;
                        C[ttt] = c20; C[ttt+1] = c21; C[ttt+2] = c22;
                    }

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

int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib, end;
    for (ib=0; ib<n; ib+=b) {
        end = ib+b-1;
        //PLU
        double max;
        int i, maxind, t, j, k;
        for (i=ib; i<=end && i<n; i++) {
            // ensure A can use LU decomposition
            maxind = i;
            max = fabs(A[i*n+i]);
            // printf("初始 max 的值为 %e", max);
            for (t=i+1; t<n; t++)
                if (fabs(A[t*n+i])>max) {
                    maxind = t;
                    max = fabs(A[t*n+i]);
                    // printf("max 的值为 %e", max);
                }
            if (max == 0) {
                printf("%s/n","LU factoration failed: coefficient matirx is singular..？");
                return -1;
            }
            else {
                if (maxind != i) {
                    int temps, index;
                    double *tempv=(double *) malloc(sizeof(double) * n);
                    temps = ipiv[i];
                    ipiv[i] = ipiv[maxind];
                    ipiv[maxind] = temps;
                    // swap rows
                    for (index = 0;index < n; index++){ // index=i
                        tempv[index] = A[i*n+index];
                        A[i*n+index] = A[maxind*n+index];
                        A[maxind*n+index] = tempv[index];
                    } 
                }
            }
            // factorization
            for (j = i+1; j<n; j++) {
                A[j*n+i] = A[j*n+i]/A[i*n+i];
                for (k=i+1; k<n && k<=end; k++)
                    A[j*n+k] -= A[j*n+i]*A[i*n+k];
            }
            // PIL done, blue&brown
        }
        //pink part: update U12
        double sum;
        for (i=ib; i<=end && i<n; i++) {
            for (j=end+1; j<n; j++) {
                sum = 0.0;
                for(k=ib; k<i; k++) {
                    sum += A[i*n+k] * A[k*n+j];
                }
                A[i*n+j] = A[i*n+j] - sum; // /1
            }
        }
 
        // // green part: use BLAS3, the most effecient one

        for (i = ib+b; i < n; i += b){
            for (j = ib+b; j < n; j += b){
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }
    return 0;
}