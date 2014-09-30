/* To Compile:
 * 'mex GEN_binomial_coeffs.c' 
 * CALL: coeffs=GEN_binomial_coeffs(N).*/

#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define max(A, B)       ((A) > (B) ? (A) : (B))
#define min(A, B)       ((A) < (B) ? (A) : (B))

/*input variables:*/
int N;

/*main function:*/
void Run(double *, int);

void Run(double *coeffs, int N) {

  int r;
  coeffs[0]=1;
  for (r=1;r<=N;r++) {
	coeffs[r]=(N+1-r)*(coeffs[r-1]/r);	
  }
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{ 
	double *coeffs;
  	N = mxGetPr(prhs[0])[0];
  	plhs[0] = mxCreateDoubleMatrix(N+1,1,mxREAL);
  	coeffs = mxGetPr(plhs[0]);
  	Run(coeffs,N);
  	return;
}
