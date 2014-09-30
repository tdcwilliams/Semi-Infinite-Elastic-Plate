/* To Compile:
 * 'mex Gaulag0.c' (on a machine with gsl installed) */
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define EPS 3.0e-14
#define MAXIT 50 

int M,N;

void Ein(double *,double *,double *);

void Ein(double *yy, double *zz, double *test) {

int i, its, j, r;
double y, z, dy;

 for (i=1;i<=M;i++) {
   for (j=1;j<=N;j++) {
     r=i-1+M*(j-1);
     z=zz[r]; dy=z; y=dy; test[0]=dy;
     for (its=2;its<=MAXIT;its++) {
	dy=-z*(its-1)*dy*1.0/(its*its);
	y=y+dy; test[its-1]=dy;
	if(fabs(dy)<EPS) break;
     }
     yy[r]=y;
   }
 }
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{
  double *yy, *zz, *test;
  zz=mxGetPr(prhs[0]);
  M=mxGetM(prhs[0]);
  N=mxGetN(prhs[0]);

  plhs[0]=mxCreateDoubleMatrix(M,N,mxREAL);
  yy = mxGetPr(plhs[0]);
  plhs[1]=mxCreateDoubleMatrix(MAXIT,1,mxREAL);
  test = mxGetPr(plhs[1]);

  Ein(yy, zz, test);

  return;
}
