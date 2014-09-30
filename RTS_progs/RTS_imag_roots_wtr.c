/* To Compile:
 * 'mex RTS_imag_roots_wtr.c' (on a machine with gsl installed) */
/*CALL: x=imag_roots_wtr(lam,H,N)
	->finds the N0-th to N1-th imaginary roots of the
	  dispersion relation for ice*/
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979
#define EPS 1.0e-8
#define MAXIT 50

/*declare input variables*/
double lam,H;
int N;

/*construct main function*/
void imag_roots_wtr(double *, double, double, int);

void imag_roots_wtr(double *x, double lam, double H, int N) {

int its, j;
double p, dp, w, w1, i0, i1;

 for (j=1;j<=N;j++) {
  /*initial guess:*/
  w=j*PI; i1=w; i0=w-PI/2.0;
  /*find a root using N-R:*/
  for (its=1;its<=MAXIT;its++) {
    p=lam*w*sin(w)+H*cos(w);
    dp=(lam-H)*sin(w)+lam*w*cos(w);
    w1=w;
    w=w1-p/dp;
    if (fabs(w-w1) <= EPS) break;
  }
  /*check root has converged & is in correct interval:*/
  if ( (its < MAXIT) && ((w-i0+EPS)*(i1+EPS-w)>=0.0) ) {
    x[j-1]=w;
  } else {
    x[j-1]=0.0;
  }
 }
  return;
}

/*define inputs and outputs*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  double *x;
  lam = *mxGetPr(prhs[0]);
  H = *mxGetPr(prhs[1]);
  N = *mxGetPr(prhs[2]);
  plhs[0]=mxCreateDoubleMatrix(N,1,mxREAL);
  x = mxGetPr(plhs[0]);

  imag_roots_wtr(x,lam,H,N);

  return;
}
