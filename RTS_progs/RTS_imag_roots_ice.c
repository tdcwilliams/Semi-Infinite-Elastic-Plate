/* To Compile:
 * 'mex RTS_imag_roots_ice.c' (on a machine with gsl installed) */
/*CALL: x=imag_roots_iceC(del,H,N0,N1)
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
double del,H;
int N0,N1;

/*construct main function*/
void imag_roots_iceC(double *, double, double, int, int);

void imag_roots_iceC(double *x, double del, double H, int N0, int N1) {

int its, j;
double p, dp, w, w1, i0, i1;
double w4, H4;

 for (j=N0;j<=N1;j++) {
  /*initial guess*/
  w=j*PI; i1=w; i0=i1-PI/2.0;
  if (j==N0) {
    i0=i1-PI;
  }
  for (its=1;its<=MAXIT;its++) {
    w4=w*w*w*w; H4=H*H*H*H;
    p=(w4+del*H4)*w*sin(w)+H4*H*cos(w);
    dp=( 5*w4+H4*(del-H) )*sin(w)+(w4+del*H4)*w*cos(w);
    w1=w;
    w=w1-p/dp;
    if (fabs(w-w1) <= EPS) break;
  }
  /*check root has converged & is in correct interval:*/
  if ( (its < MAXIT) && ((w-i0+EPS)*(i1+EPS-w)>=0.0) ) {
    x[j-N0]=w;
  } else {
    x[j-N0]=0.0;
  }
 }
  return;
}


 
/*define inputs and outputs*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{
  double *x;
  del = *mxGetPr(prhs[0]);
  H = *mxGetPr(prhs[1]);
  N0 = *mxGetPr(prhs[2]);
  N1 = *mxGetPr(prhs[3]);
  plhs[0]=mxCreateDoubleMatrix(N1+1-N0,1,mxREAL);
  x = mxGetPr(plhs[0]);

  imag_roots_iceC(x,del,H,N0,N1);

  return;
}
