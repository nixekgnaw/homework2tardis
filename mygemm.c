#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

void dgemm0(const double *A, const double *B, double *C, const int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                C[i * n + j] += A[i * n + k] * B[k * n + j];
}

void dgemm1(const double *A, const double *B, double *C, const int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            register double r = C[i * n + j];
            for (int k = 0; k < n; k++)
                r += A[i * n + k] * B[k * n + j];
            C[i * n + j] = r;
        }
}

void dgemm2(const double *A, const double *B, double *C, const int n)
{
    int i, j, k;
    for (i = 0; i < n; i += 2)
        for (j = 0; j < n; j += 2)
        {
            register int t = i * n + j;
            register int tt = t + n;
            register double c00 = C[t];
            register double c01 = C[t + 1];
            register double c10 = C[tt];
            register double c11 = C[tt + 1];
            for (k = 0; k < n; k += 2)
            { /* 2 by 2 mini matrix multiplication using registers*/
                register int ta = i * n + k;
                register int tta = ta + n;
                register int tb = k * n + j;
                register int ttb = tb + n;
                register double a00 = A[ta];
                register double a01 = A[ta + 1];
                register double a10 = A[tta];
                register double a11 = A[tta + 1];
                register double b00 = B[tb];
                register double b01 = B[tb + 1];
                register double b10 = B[ttb];
                register double b11 = B[ttb + 1];
                c00 += a00 * b00 + a01 * b10;
                c01 += a00 * b01 + a01 * b11;
                c10 += a10 * b00 + a11 * b10;
                c11 += a10 * b01 + a11 * b11;
            }
            C[t] = c00;
            C[t + 1] = c01;
            C[tt] = c10;
            C[tt + 1] = c11;
        }
}

void dgemm3(const double *A, const double *B, double *C, const int n)
{
    for (int i = 0; i < n; i += 3)
    {
        for (int j = 0; j < n; j += 3)
        {
            register int t = i * n + j;
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
            for (int k = 0; k < n; k += 3)
            {
                register int ta = i * n + k;
                register int tta = ta + n;
                register int ttta = tta + n;
                register int tb = k * n + j;
                register int ttb = tb + n;
                register int tttb = ttb + n;
                /*
                * C00 = A00*B00 + A01*B10 + A02*B20
                * C01 = A00*B01 + A01*B11 + A02*B21
                * C02 = A00*B02 + A01*B12 + A02*B22
                *
                * C10 = A10*B00 + A11*B10 + A12*B20;
                * C11 = A10*B01 + A11*B11 + A12*B21;
                * C12 = A10*B02 + A11*B12 + A12*B22;
                *
                * C20 = A20*B00 + A21*B10 + A22*B20;
                * C21 = A20*B01 + A21*B11 + A22*B21;
                * C22 = A20*B02 + A21*B12 + A22*B22;
                */
                register double r1 = A[ta];     // a00
                register double r2 = A[tta];    // a10
                register double r3 = A[ttta];   // a20
                register double r4 = B[tb];     // b00
                register double r5 = B[tb + 1]; // b01
                register double r6 = B[tb + 2]; // b02
                c00 += r1 * r4;
                c01 += r1 * r5;
                c02 += r1 * r6;
                c10 += r2 * r4;
                c11 += r2 * r5;
                c12 += r2 * r6;
                c20 += r3 * r4;
                c21 += r3 * r5;
                c22 += r3 * r6;
                r1 = A[ta + 1];
                r2 = A[tta + 1];
                r3 = A[ttta + 1];
                r4 = B[ttb];
                r5 = B[ttb + 1];
                r6 = B[ttb + 2];
                c00 += r1 * r4;
                c01 += r1 * r5;
                c02 += r1 * r6;
                c10 += r2 * r4;
                c11 += r2 * r5;
                c12 += r2 * r6;
                c20 += r3 * r4;
                c21 += r3 * r5;
                c22 += r3 * r6;
                r1 = A[ta + 2];
                r2 = A[tta + 2];
                r3 = A[ttta + 2];
                r4 = B[tttb];
                r5 = B[tttb + 1];
                r6 = B[tttb + 2];
                c00 += r1 * r4;
                c01 += r1 * r5;
                c02 += r1 * r6;
                c10 += r2 * r4;
                c11 += r2 * r5;
                c12 += r2 * r6;
                c20 += r3 * r4;
                c21 += r3 * r5;
                c22 += r3 * r6;
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

void ijk(const double *A, const double *B, double *C, const int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            register double sum = 0.0;
            for (int k = 0; k < n; k++)
                sum += A[i * n + k] * B[k * n + j];
            C[i * n + j] = sum;
        }
    }
}

void bijk(const double *A, const double *B, double *C, const int n, const int b)
{
    for (int i = 0; i < n; i += b)
        for (int j = 0; j < n; j += b)
            for (int k = 0; k < n; k += b)
                for (int i1 = i; i1 < i + b; i1++)
                {
                    for (int j1 = j; j1 < j + b; j1++)
                    {
                        register double sum = 0.0;
                        for (int k1 = k; k1 < k + b; k1++)
                            sum += A[i1 * n + k1] * B[k1 * n + j1];
                        C[i1 * n + j1] = sum;
                    }
                }
}

void jik(const double *A, const double *B, double *C, const int n)
{
    /* jik */
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            register double sum = 0.0;
            for (int k = 0; k < n; k++)
                sum += A[i * n + k] * B[k * n + j];
            C[i * n + j] = sum;
        }
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b)
{
    for (int j = 0; j < n; j += b)
        for (int i = 0; i < n; i += b)
            for (int k = 0; k < n; k += b)
                for (int j1 = j; j1 < j + b; j1++)
                {
                    for (int i1 = i; i1 < i + b; i1++)
                    {
                        register double sum = 0.0;
                        for (int k1 = k; k1 < k + b; k1++)
                            sum += A[i1 * n + k1] * B[k1 * n + j1];
                        C[i1 * n + j1] = sum;
                    }
                }
}

void kij(const double *A, const double *B, double *C, const int n)
{
    /* kij */
    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < n; i++)
        {
            register double r = A[i * n + k];
            for (int j = 0; j < n; j++)
                C[i * n + j] += r * B[k * n + j];
        }
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b)
{
    for (int k = 0; k < n; k += b)
        for (int i = 0; i < n; i += b)
            for (int j = 0; j < n; j += b)
                for (int k1 = k; k1 < k + b; k1++)
                {
                    for (int i1 = i; i1 < i + b; i1++)
                    {
                        register double r = A[i1 * n + k1];
                        for (int j1 = j; j1 < j + b; j1++)
                            C[i1 * n + j1] += r * B[k1 * n + j1];
                    }
                }
}

void ikj(const double *A, const double *B, double *C, const int n)
{
    /* ikj */
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < n; k++)
        {
            register double r = A[i * n + k];
            for (int j = 0; j < n; j++)
                C[i * n + j] += r * B[k * n + j];
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b)
{
    for (int i = 0; i < n; i += b)
        for (int k = 0; k < n; k += b)
            for (int j = 0; j < n; j += b)
                for (int i1 = i; i1 < i + b; i1++)
                {
                    for (int k1 = k; k1 < k + b; k1++)
                    {
                        register double r = A[i1 * n + k1];
                        for (int j1 = j; j1 < j + b; j1++)
                            C[i1 * n + j1] += r * B[k1 * n + j1];
                    }
                }
}

void jki(const double *A, const double *B, double *C, const int n)
{
    /* jki */
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < n; k++)
        {
            register double r = B[k * n + j];
            for (int i = 0; i < n; i++)
                C[i * n + j] += A[i * n + k] * r;
        }
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b)
{
    for (int j = 0; j < n; j += b)
        for (int k = 0; k < n; k += b)
            for (int i = 0; i < n; i += b)
                for (int j1 = j; j1 < j + b; j1++)
                {
                    for (int k1 = k; k1 < k + b; k1++)
                    {
                        register double r = B[k1 * n + j1];
                        for (int i1 = i; i1 < i + b; i1++)
                            C[i1 * n + j1] += A[i1 * n + k1] * r;
                    }
                }
}

void kji(const double *A, const double *B, double *C, const int n)
{
    /* kji */
    for (int k = 0; k < n; k++)
    {
        for (int j = 0; j < n; j++)
        {
            register double r = B[k * n + j];
            for (int i = 0; i < n; i++)
                C[i * n + j] += A[i * n + k] * r;
        }
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b)
{
    for (int k = 0; k < n; k += b)
        for (int j = 0; j < n; j += b)
            for (int i = 0; i < n; i += b)
                for (int k1 = k; k1 < k + b; k1++)
                {
                    for (int j1 = j; j1 < j + b; j1++)
                    {
                        register double r = B[k1 * n + j1];
                        for (int i1 = i; i1 < i + b; i1++)
                            C[i1 * n + j1] += A[i1 * n + k1] * r;
                    }
                }
}

void optimal(const double *A, const double *B, double *C, const int n, const int b)
{
    
}