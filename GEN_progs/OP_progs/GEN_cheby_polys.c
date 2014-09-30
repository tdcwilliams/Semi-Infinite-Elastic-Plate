/* To Compile:
 * 'mex GEN_cheby_polys.c' */
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define max(A, B)       ((A) > (B) ? (A) : (B))
#define min(A, B)       ((A) < (B) ? (A) : (B))
#define EPS 

int N;

void fun(double *, double *, int);

void fun(double *T, double *x, int n_cheby) {
  int i,j,r;
  for(i=1;i<=N;i++){
	T[i-1]=1; T[i-1+N]=x[i-1];
     for(j=2;j<=n_cheby;j++){
	T[i-1+j*N]=2*x[i-1]*T[i-1+(j-1)*N]-T[i-1+(j-2)*N];
     }
  }
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
/*CALL: sum=GEN_run_sum(vec), where vec is a REAL col vector of length N.
	Function returns an (N+1)-vector containing running sums of 
	the elements of vec.*/
{
  double *T, *x, *P;
  int n_cheby;
/* int N, n_cheby;*/

  x = mxGetPr(prhs[0]);
  N = mxGetM(prhs[0]);/*=length(x)*/

  n_cheby=mxGetPr(prhs[1])[0];/*T_{n_cheby} is the highest degree chebyshev
				polynomial required.*/
  plhs[0] = mxCreateDoubleMatrix(N,n_cheby+1,mxREAL);
  T = mxGetPr(plhs[0]);

  fun(T,x,n_cheby);
  return;
}
