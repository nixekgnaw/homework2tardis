#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/
void calLastRowandColumn(const double* A, const double* B, double* C, const int n);

void dgemm0(const double* A, const double* B, double* C, const int n)
{
	int i,j,k;
	for(i =0;i<n;i++){
		for (j = 0; j < n ; j++){
			for(k=0;k<n;k++){
				C[i*n+j] += A[i*n+k] * B[k*n+j];
			}
		}
	}
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
	int i,j,k;
	for (i=0; i<n; i++){
		for (j=0; j<n; j++) {
			register double r = C[i*n+j] ;
			for (k=0; k<n; k++){
				r += A[i*n+k] * B[k*n+j];
			}
			C[i*n+j] = r;
		}
	}
}

void dgemm2(const double *A, const double *B, double *C, const int n) 
{
	int i,j,k;
	for (i = 0; i < n; i+=2) {
		for (j = 0; j < n; j += 2) {
			register int t = i*n+j;
			register int tt = t+n;
			register double c00 = C[t];
			register double c01 = C[t+1];
			register double c10 = C[tt];
			register double c11 = C[tt+1];

			for (k = 0; k < n; k += 2) {
				/* 2 by 2 mini matrix multiplication using registers*/
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
			C[t+1] = c01;
			C[tt] = c10;
			C[tt+1] = c11;
		}
	}
}

void dgemm3(const double *A, const double *B, double *C, const int n) 
{

	int i,j,k;
	int boundary = n-(n%3);
	for (i = 0; i <boundary;i+=3) {
		for(j=0; j<boundary;j+=3){
			register int t = i*n+j; //row 0
			register int tt = t+n;  //row 1
			register int ttt = tt+n; //row 2
			//define 9 register for c result.
			register double c00 = C[t];
			register double c01 = C[t+1];
			register double c02 = C[t+2];

			register double c10 = C[tt];
			register double c11 = C[tt+1];
			register double c12 = C[tt+2];

			register double c20 = C[ttt];
			register double c21 = C[ttt+1];
			register double c22 = C[ttt+2];

			register int ta ;
			register int tta ;
			register int ttta ;

			register int tb;
			register int ttb;
			register int tttb;
			//define 6 register for a and b.
			register double a00;
			register double a10;
			register double a20;

			register double b00;
			register double b01;
			register double b02;

			for(k = 0; k<boundary; k+=3){
				ta = i*n+k;  // a row 0
				tta = ta+n;  // a row 1
				ttta = tta+n; // a row 2

				tb = k*n+j;   // b column 0
				ttb = tb +n;   // b column 1
				tttb = ttb +n;   // b column 2

				//start: calculate value of one block in one inner loop
				a00 = A[ta];  //row1
				a10 = A[tta]; // row2
				a20 = A[ttta];//row3

				b00 = B[tb];  // column 1
				b01 = B[tb+1];// column 2
				b02 = B[tb+2];//column 3

				//first
				c00 += a00 * b00;
				c01 += a00 * b01;
				c02 += a00 * b02;

				c10 += a10 * b00;
				c11 += a10 * b01;
				c12 += a10 * b02;

				c20 += a20 * b00;
				c21 += a20 * b01;
				c22 += a20 * b02;

				//second
				a00 = A[ta+1];
				a10 = A[tta+1];
				a20 = A[ttta+1];

				b00 = B[ttb];
				b01 = B[ttb+1];
				b02 = B[ttb+2];

				c00 += a00 * b00;
				c01 += a00 * b01;
				c02 += a00 * b02;

				c10 += a10 * b00;
				c11 += a10 * b01;
				c12 += a10 * b02;

				c20 += a20 * b00;
				c21 += a20 * b01;
				c22 += a20 * b02;

				//third
				a00 = A[ta+2];
				a10 = A[tta+2];
				a20 = A[ttta+2];

				b00 = B[tttb];
				b01 = B[tttb+1];
				b02 = B[tttb+2];

				c00 += a00 * b00;
				c01 += a00 * b01;
				c02 += a00 * b02;

				c10 += a10 * b00;
				c11 += a10 * b01;
				c12 += a10 * b02;

				c20 += a20 * b00;
				c21 += a20 * b01;
				c22 += a20 * b02;

			}
			//start to calculate the vaule from boundary to n
			for(;k<n;k++){
				ta = i*n+k;  // a row 0
				tta = ta+n;  // a row 1
				ttta = tta+n; // a row 2

				tb = k*n+j;   // b column 0

				a00 = A[ta];  //row1
				a10 = A[tta]; // row2
				a20 = A[ttta];//row3

				b00 = B[tb];  // column 1
				b01 = B[tb+1];// column 2
				b02 = B[tb+2];//column 3

				c00 += a00 * b00;
				c01 += a00 * b01;
				c02 += a00 * b02;

				c10 += a10 * b00;
				c11 += a10 * b01;
				c12 += a10 * b02;

				c20 += a20 * b00;
				c21 += a20 * b01;
				c22 += a20 * b02;
			}
			//end

			C[t] = c00;
			C[t+1] = c01;
			C[t+2] = c02;

			C[tt] = c10;
			C[tt+1] = c11;
			C[tt+2] = c12;

			C[ttt] = c20;
			C[ttt+1] = c21;
			C[ttt+2] = c22;
			//end 
		}  
	}
	//cal last row and column
	calLastRowandColumn(A,B,C,n);
	

}

void calLastRowandColumn(const double* A, const double* B, double* C, const int n){
	int remainAmount = n % 3;
	if(remainAmount > 0){
		int i,j,k;
		register double result;
		for(i = n - remainAmount;i<n;i++) {
			for (j = 0; j < n; j++) {
				result = C[i*n+j];
				for (k = 0; k < n; k++) {
					result += A[i * n + k] * B[k * n + j];
				}
				C[i * n + j] = result;
			}
		}

		for(i = 0;i<n-remainAmount;i++){
			for(j = n-remainAmount; j<n;j++){
				result = C[i*n+j];
				for(k=0;k<n;k++){
					result += A[i * n + k] * B[k * n + j];
				}
				C[i * n + j] = result;
			}
		}
	}
}

void ijk(const double *A, const double *B, double *C, const int n) 
{
	int i,j,k;
	for (i=0; i<n; i++)  {
		for (j=0; j<n; j++) {
			register double sum = C[i*n+j];
			for (k=0; k<n; k++) 
				sum += A[i*n+k] * B[k*n+j];
			C[i*n+j] = sum;
		}
	} 

}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
	/* ijk – blocked version algorithm*/ 
	int i,j,k,i1,j1,k1;
	for (i = 0; i < n; i+=b) 
		for (j = 0; j < n; j+=b) 
			for (k = 0; k < n; k+=b)
				/* B x B mini matrix multiplications */ 
				for (i1 = i; i1 < i+b && i1<n; i1++)
					for (j1 = j; j1 < j+b && j1<n; j1++) { 
						register double r=C[i1*n+j1]; 
							for (k1 = k; k1 < k+b && k1<n; k1++) 
								r += A[i1*n + k1]*B[k1*n + j1]; 
						C[i1*n+j1]=r; 
					}
} 


void jik(const double *A, const double *B, double *C, const int n) 
{
	int i,j,k;
	for (j=0; j<n; j++)  {
		for (i=0; i<n; i++) {
			register double sum = C[i*n+j];
			for (k=0; k<n; k++) 
				sum += A[i*n+k] * B[k*n+j];
			C[i*n+j] = sum;
		}
	} 

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i,j,k,i1,j1,k1;
	for (j = 0; j < n; j+=b) 
		for (i = 0; i < n; i+=b) 
			for (k = 0; k < n; k+=b)
				/* B x B mini matrix multiplications */ 
				for (j1 = j; j1 < j+b && j1<n; j1++)
					for (i1 = i; i1 < i+b && i1<n; i1++) { 
						register double r=C[i1*n+j1]; 
							for (k1 = k; k1 < k+b && k1<n; k1++) 
								r += A[i1*n + k1]*B[k1*n + j1]; 
						C[i1*n+j1]=r; 
					}
					
}

void kij(const double *A, const double *B, double *C, const int n) 
{

	int i,j,k;
	for (k=0; k<n; k++) {
		for (i=0; i<n; i++) {
			register double r = A[i*n+k];
			for (j=0; j<n; j++)
				C[i*n+j] += r * B[k*n+j];
		}
	}
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i,j,k,i1,j1,k1;
	 for (k = 0; k < n; k+=b)
		for (i = 0; i < n; i+=b) 
			for (j = 0; j < n; j+=b)
				/* B x B mini matrix multiplications */ 
				 for (k1 = k; k1 < k+b && k1<n; k1++)
					for (i1 = i; i1 < i+b && i1<n; i1++) { 
						register double r=A[i1*n+k1]; 
							for (j1 = j; j1 < j+b && j1<n; j1++)
								C[i1*n+j1] += r*B[k1*n + j1]; 
					}
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
	int i,j,k;
	for (i=0; i<n; i++) {
		for (k=0; k<n; k++) {
			register double  r = A[i*n+k];
			for (j=0; j<n; j++)
				C[i*n+j] += r * B[k*n+j];
		}
	}

}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i,j,k,i1,j1,k1;
	  for (i = 0; i < n; i+=b) 
		for (k = 0; k < n; k+=b) 
			for (j = 0; j < n; j+=b)
				/* B x B mini matrix multiplications */ 
				  for (i1 = i; i1 < i+b && i1<n; i1++)
					for (k1 = k; k1 < k+b && k1<n; k1++) { 
						register double r=A[i1*n+k1]; 
							for (j1 = j; j1 < j+b && j1<n; j1++)
								C[i1*n+j1] += r*B[k1*n + j1]; 
					}
}

void jki(const double *A, const double *B, double *C, const int n) 
{
	int i,j,k;
	for (j=0; j<n; j++) {
		for (k=0; k<n; k++) {
			register double r = B[k*n+j];
			for (i=0; i<n; i++)
				C[i*n+j] += A[i*n+k] * r;
		}
	}
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i,j,k,i1,j1,k1;
	   for (j = 0; j < n; j+=b)
		for (k = 0; k < n; k+=b) 
			for (i = 0; i < n; i+=b)
				/* B x B mini matrix multiplications */ 
				   for (j1 = j; j1 < j+b && j1<n; j1++)
					for (k1 = k; k1 < k+b && k1<n; k1++) { 
						register double r=B[k1*n+j1]; 
							for (i1 = i; i1 < i+b && i1<n; i1++)
								C[i1*n+j1] += A[i1*n+k1]*r; 
					}

}

void kji(const double *A, const double *B, double *C, const int n) 
{
	int i,j,k;
	for (k=0; k<n; k++) {
		for (j=0; j<n; j++) {
			register double	r = B[k*n+j];
			for (i=0; i<n; i++)
				C[i*n+j] += A[i*n+k] * r;
		}
	}
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
	int i,j,k,i1,j1,k1;
	    for (k = 0; k < n; k+=b)
		 for (j = 0; j < n; j+=b)
			for (i = 0; i < n; i+=b)
				/* B x B mini matrix multiplications */ 
				    for (k1 = k; k1 < k+b && k1<n; k1++)
					 for (j1 = j; j1 < j+b && j1<n; j1++){ 
						register double r=B[k1*n+j1]; 
							for (i1 = i; i1 < i+b && i1<n; i1++)
								C[i1*n+j1] += A[i1*n+k1]*r; 
					 }
}

void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
	int i,j,k,i1,j1,k1;
	for (i = 0; i < n; i+=b) 
		for (j = 0; j < n; j+=b) 
			for (k = 0; k < n; k+=b)
				/* B x B mini matrix multiplications */ 
				for (i1 = i; i1 < i+b &&i1<n; i1+=3)
					for (j1 = j; j1<j+b && j1<n; j1+=3) { 
			register int t = i1*n+j1; //row 0
			register int tt = t+n;  //row 1
			register int ttt = tt+n; //row 2
			//define 9 register for c result.
			register double c00 = C[t];
			register double c01 = C[t+1];
			register double c02 = C[t+2];

			register double c10 = C[tt];
			register double c11 = C[tt+1];
			register double c12 = C[tt+2];

			register double c20 = C[ttt];
			register double c21 = C[ttt+1];
			register double c22 = C[ttt+2];

			register int ta ;
			register int tta ;
			register int ttta ;

			register int tb;
			register int ttb;
			register int tttb;
			//define 6 register for a and b.
			register double a00;
			register double a10;
			register double a20;

			register double b00;
			register double b01;
			register double b02;
			for (k1 = k; k1 < k+b && k1<n; k1+=3){
				ta = i1*n+k1;  // a row 0
				tta = ta+n;  // a row 1
				ttta = tta+n; // a row 2

				tb = k1*n+j1;   // b column 0
				ttb = tb +n;   // b column 1
				tttb = ttb +n;   // b column 2

				//start: calculate value of one block in one inner loop
				a00 = A[ta];  //row1
				a10 = A[tta]; // row2
				a20 = A[ttta];//row3

				b00 = B[tb];  // column 1
				b01 = B[tb+1];// column 2
				b02 = B[tb+2];//column 3

				//first
				c00 += a00 * b00;
				c01 += a00 * b01;
				c02 += a00 * b02;

				c10 += a10 * b00;
				c11 += a10 * b01;
				c12 += a10 * b02;

				c20 += a20 * b00;
				c21 += a20 * b01;
				c22 += a20 * b02;

				//second
				a00 = A[ta+1];
				a10 = A[tta+1];
				a20 = A[ttta+1];

				b00 = B[ttb];
				b01 = B[ttb+1];
				b02 = B[ttb+2];

				c00 += a00 * b00;
				c01 += a00 * b01;
				c02 += a00 * b02;

				c10 += a10 * b00;
				c11 += a10 * b01;
				c12 += a10 * b02;

				c20 += a20 * b00;
				c21 += a20 * b01;
				c22 += a20 * b02;

				//third
				a00 = A[ta+2];
				a10 = A[tta+2];
				a20 = A[ttta+2];

				b00 = B[tttb];
				b01 = B[tttb+1];
				b02 = B[tttb+2];

				c00 += a00 * b00;
				c01 += a00 * b01;
				c02 += a00 * b02;

				c10 += a10 * b00;
				c11 += a10 * b01;
				c12 += a10 * b02;

				c20 += a20 * b00;
				c21 += a20 * b01;
				c22 += a20 * b02;

			}
			C[t] = c00;
			C[t+1] = c01;
			C[t+2] = c02;

			C[tt] = c10;
			C[tt+1] = c11;
			C[tt+2] = c12;

			C[ttt] = c20;
			C[ttt+1] = c21;
			C[ttt+2] = c22;
		}
}
