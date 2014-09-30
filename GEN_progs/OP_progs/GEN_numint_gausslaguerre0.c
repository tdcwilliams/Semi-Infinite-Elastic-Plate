/* To Compile:
 * 'mex Gaulag0.c' (on a machine with gsl installed) */
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>

#define EPS 3.0e-14
#define MAXIT 10 

int n;
double alf;

void Gaulag0(double *,double *,double , int );

void Gaulag0(double *x, double *w, double alf, int n) {

int i, its, j;
float ai;
double p1, p2, p3, pp, z, z0, z1;

 for (i=1;i<=n;i++) {
  if (i==1) {
    z = (1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
  } else if (i==2) {
    z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
  } else {
    ai = i-2;
    z0 = ((1.0+2.55*ai)/(1.9*ai) + 1.26*ai*alf/(1.0+3.5*ai));
    z += z0*(z-x[i-3])/(1.0+.3*alf);
  }
  for (its=1;its<=MAXIT;its++) {
    p1 = 1.0;
    p2 = 0.0;
    for (j=1;j<=n;j++) {
      p3 = p2;
      p2 = p1;
      p1 = ( (2*j-1+alf-z)*p2 - (j-1+alf)*p3 )/j;
    }
    pp=(n*p1-(n+alf)*p2)/z;
    z1=z;
    z=z1-p1/pp;
    if (fabs(z-z1) <= EPS) break;
  }
  x[i-1]=z;
  w[i-1]=1/pp/p2;
 }
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{
  double *x, *w;
  alf = *mxGetPr(prhs[0]);
  n = *mxGetPr(prhs[1]);
  plhs[0]=mxCreateDoubleMatrix(n,1,mxREAL);
  plhs[1]=mxCreateDoubleMatrix(n,1,mxREAL);
  x = mxGetPr(plhs[0]);
  w = mxGetPr(plhs[1]);

  Gaulag0(x ,w, alf, n);

  return;
}
