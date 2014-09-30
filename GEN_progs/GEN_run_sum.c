/* To Compile:
 * 'mex GEN_run_sum.c' */
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define max(A, B)       ((A) > (B) ? (A) : (B))
#define min(A, B)       ((A) < (B) ? (A) : (B))
#define EPS 

int N;

void fun(double *, double *);

void fun(double *sum, double *vec) {
  int j;
  sum[0]=0;
  for(j=1;j<=N;j++){
	sum[j]=sum[j-1]+vec[j-1];
  }
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
/*CALL: sum=GEN_run_sum(vec), where vec is a REAL col vector of length N.
	Function returns an (N+1)-vector containing running sums of 
	the elements of vec.*/
{
  double *sum, *vec;
  vec = mxGetPr(prhs[0]);
  N = mxGetM(prhs[0]);/*=length(vec)*/
  plhs[0] = mxCreateDoubleMatrix(N+1,1,mxREAL);
  sum = mxGetPr(plhs[0]);

  fun(sum,vec);
  return;
}
