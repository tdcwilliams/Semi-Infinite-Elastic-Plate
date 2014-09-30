/* To Compile:
 * 'mex GEN_numdiff.c' */
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define max(A, B)       ((A) > (B) ? (A) : (B))
#define min(A, B)       ((A) < (B) ? (A) : (B))
#define EPS 

void run_fun(double *, int);

void run_fun(double *d_op, int N) {
  int i, j, r, s, t, sgn;
  int j_m, j_p, r_m, r_p;
  /*N=N_pr[0];*/
  /*do 1st & last rows 1st:*/
  for (s=1;s<=2;s++){
    for(t=1;t<=2;t++){
      i=1+(s-1)*N;/*s==1 -> i=1; s==2 -> i=N+1;*/
      j=t+(s-1)*(N-1);/*s==1 -> j=1,2; s==2 -> j=N,N+1;*/
      sgn=1-2*(t==1);/*t==1 -> -1,t==2 -> +1;*/
      r=i-1+(j-1)*(N+1); d_op[r]=2*sgn;
    }
  }
  /*now do remaining rows:*/
  for(i=2;i<=N;i++){
    j_m=i-1; j_p=i+1;
    r_m=i-1+(j_m-1)*(N+1); r_p=i-1+(j_p-1)*(N+1);
    d_op[r_m]=-1; d_op[r_p]=1;
  }
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
/*CALL: d_op=GEN_numdiff(N), where N is an integer.
	Function returns an (N+1)x(N+1) matrix, d_op; 
	given an (N+1)-vector f, df=d_op*f/2/Del.*/
{
  double *N_pr, *d_op;
  int N;
  N_pr = mxGetPr(prhs[0]);
  N=N_pr[0];/*N=N_pr[0]=no of panels*/
  plhs[0] = mxCreateDoubleMatrix(N+1,N+1,mxREAL);
  d_op = mxGetPr(plhs[0]);

  run_fun(d_op,N);
  return;
}
