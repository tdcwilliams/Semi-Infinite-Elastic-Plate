/* To Compile:
 * 'mex GEN_DAD.c' 
 * CALL: M=GEN_DAD(v1,A,v2); 
 * M=diag(v1)*A*diag(v2), i.e. M(i,j)=v1(i)*A(i,j)*v2(j).*/

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
double *v1_re, *v1_im, *A_re, *A_im, *v2_re, *v2_im;

/*main function:*/
void Run(double *,double *,double *,double *,double *,double *,double *,double *,double *);

void Run(double *M_re,double *M_im, double *v1_re,double *v1_im, double *A_re,double *A_im, double *v2_re,double *v2_im, double *test) {

  int i, j, r;
  double a, b, c, d, e, f;
  for (i=1;i<=M;i++) {
    for (j=1;j<=N;j++) {
	r=(j-1)*M+i-1;
	a=v1_re[i-1]; b=v1_im[i-1];
	c=A_re[r]; d=A_im[r];
	e=v2_re[j-1]; f=v2_im[j-1];
	M_re[r]=(a*c-b*d)*e-(b*c+a*d)*f;
	M_im[r]=(a*c-b*d)*f+(b*c+a*d)*e;
    }
  }
	i=1;j=1; r=(j-1)*M+i-1;
 a=10; b=11; c=1; d=2; e=3; f=4;
	a=v1_re[i-1]; b=v1_im[i-1];/*ok*/
	c=A_re[r]; d=A_im[r];/*ok*/
	e=v2_re[j-1]; f=v2_im[j-1];/*ok*/
	test[0]=a; test[1]=b; test[2]=c; test[3]=d; test[4]=e; test[5]=f;
  test[6]=a*c-b*d;test[7]=b*c+a*d; test[8]=(a*c-b*d)*e-f*(b*c+a*d);
  test[9]=(a*c-b*d)*f-e*(b*c+a*d); test[10]=r;
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{ 
  	double *M_re, *M_im;
	double *test;
  	M = mxGetM(prhs[0]);/*length of v1*/
	N = mxGetM(prhs[2]);/*length of v2*/
	v1_re = mxGetPr(prhs[0]);
	v1_im = mxGetPi(prhs[0]);/*real and imaginary parts of v1*/
	A_re = mxGetPr(prhs[1]);
	A_im = mxGetPi(prhs[1]);/*real and imaginary parts of A*/
	v2_re = mxGetPr(prhs[2]);
	v2_im = mxGetPi(prhs[2]);/*real and imaginary parts of v2*/

  	plhs[0] = mxCreateDoubleMatrix(M,N,mxCOMPLEX);
  	M_re = mxGetPr(plhs[0]);
	M_im = mxGetPi(plhs[0]);/*real and imaginary parts of Ap11*/

	plhs[1]=mxCreateDoubleMatrix(15,1,mxREAL);
	test = mxGetPr(plhs[1]);

  	Run(M_re,M_im,v1_re,v1_im,A_re,A_im,v2_re,v2_im,test);
  	return;
}
