/* To Compile:
 * 'mex GEN_ADr.c' 
 * CALL: M=GEN_ADr(A,v) (v real);
 * M=A*diag(v), i.e. M(i,j)=A(i,j)*v(j).*/

#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
/*#include <gsl/gsl_math.h>*/
/*#include <complex_math.h>*/

#define max(A, B)       ((A) > (B) ? (A) : (B))
#define min(A, B)       ((A) < (B) ? (A) : (B))

/*input variables:*/
int M, N;
double *A_re, *A_im, *v;

/*main function:*/
void Run(double *,double *,double *,double *,double *,double *);

void Run(double *M_re,double *M_im, double *A_re,double *A_im, double *v, double *test) {

  int i, j, r;
  double a, b, c, d, e, f;
  for (i=1;i<=M;i++) {
    for (j=1;j<=N;j++) {
	r=(j-1)*M+i-1;
	M_re[r]=A_re[r]*v[j-1];
	M_im[r]=A_im[r]*v[j-1];
    }
  }
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{ 
  	double *M_re, *M_im;
	double *test;
  	M = mxGetM(prhs[0]);
	N = mxGetN(prhs[0]);/*A is MxN*/
	A_re = mxGetPr(prhs[0]);
	A_im = mxGetPi(prhs[0]);/*real and imaginary parts of A*/
	v = mxGetPr(prhs[1]);/*v is real*/

  	plhs[0] = mxCreateDoubleMatrix(M,N,mxCOMPLEX);
  	M_re = mxGetPr(plhs[0]);
	M_im = mxGetPi(plhs[0]);/*real and imaginary parts of Ap11*/

	plhs[1]=mxCreateDoubleMatrix(15,1,mxREAL);
	test = mxGetPr(plhs[1]);

  	Run(M_re,M_im,A_re,A_im,v,test);
  	return;
}
