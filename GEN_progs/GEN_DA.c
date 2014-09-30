/* To Compile:
 * 'mex GEN_DA.c' 
 * CALL: M=GEN_DA(v1,A); 
 * M=diag(v1)*A, i.e. M(i,j)=v1(i)*A(i,j).*/

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
double *v1_re, *v1_im, *A_re, *A_im;

/*main function:*/
void Run(double *,double *,double *,double *,double *,double *,double *);

void Run(double *M_re,double *M_im, double *v1_re,double *v1_im, double *A_re,double *A_im, double *test) {

  int i, j, r;
  double a, b, c, d;
  for (i=1;i<=M;i++) {
    for (j=1;j<=N;j++) {
	r=(j-1)*M+i-1;
	a=v1_re[i-1];
	b=v1_im[i-1];
	c=A_re[r];
	d=A_im[r];
	M_re[r]=a*c-b*d;
	M_im[r]=b*c+a*d;
    }
  }
/*	i=1;j=1; r=(j-1)*M+i-1;
 a=10; b=11; c=1; d=2; e=3; f=4;
	a=v1_re[i-1]; b=v1_im[i-1];*//*ok*/
/*	c=A_re[r]; d=A_im[r];*//*ok*/
/*	e=v2_re[j-1]; f=v2_im[j-1];*//*ok*/
/*	test[0]=a; test[1]=b; test[2]=c; test[3]=d; test[4]=e; test[5]=f;
  test[6]=a*c-b*d;test[7]=b*c+a*d; test[8]=(a*c-b*d)*e-f*(b*c+a*d);
  test[9]=(a*c-b*d)*f-e*(b*c+a*d); test[10]=r;*/
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{ 
  double *M_re, *M_im;
  double *test;
  M = mxGetM(prhs[1]);/*no of rows in A*/
  N = mxGetN(prhs[1]);/*no of cols in A*/
  v1_re = mxGetPr(prhs[0]);
  v1_im = mxGetPi(prhs[0]);/*real and imaginary parts of v1*/
  A_re = mxGetPr(prhs[1]);
  A_im = mxGetPi(prhs[1]);/*real and imaginary parts of A*/

  plhs[0] = mxCreateDoubleMatrix(M,N,mxCOMPLEX);
  M_re = mxGetPr(plhs[0]);
  M_im = mxGetPi(plhs[0]);/*real and imaginary parts of Ap11*/

  plhs[1]=mxCreateDoubleMatrix(15,1,mxREAL);
  test = mxGetPr(plhs[1]);

  Run(M_re,M_im,v1_re,v1_im,A_re,A_im,test);
  return;
}
