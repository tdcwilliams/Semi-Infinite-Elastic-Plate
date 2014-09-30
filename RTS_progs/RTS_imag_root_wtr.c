/* To Compile:
 * 'mex imag_root_wtr.c' (on a machine with gsl installed) */
/*CALL: x=imag_root_wtr(lam,H,i0,i1)
	->finds the imaginary root of the dispersion relation for water,
	  f=lam*w*sin(w)+H*cos(w), inside the interval [i0,i1]*/
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979
#define EPS 1.0e-8
#define MAXIT 100

/*declare input variables*/
double lam, H, i0, i1;

/*construct main function*/
void imag_root_wtr(double *, double, double, double, double);

void imag_root_wtr(double *w, double lam, double H, double i0, double i1) {

int its, j;
double p, dp, x, x0;
double xl, fl, xh, fh;
double dxold, dx;

  fl=lam*i0*sin(i0)+H*cos(i0); fh=lam*i1*sin(i1)+H*cos(i1);
if (1) {/*necessary? may save time if either endpoint is "close enough"*/
  /*check there is a root in interval*/
  if ( (fl > 0.0 && fh > 0.0) || ( fl < 0.0) && fh < 0.0 ) {
	w[0]=0.0; return;
  } 
  /*if either endpoint is "close enough", just use that for root*/
  if (fl == 0.0) {w[0] = i0; return;}
  if (fh == 0.0) {w[0] = i1; return;}
}
  /*define search direction to be in dirn of increasing p*/
  if (fl < 0.0) {
    xl=i0; xh=i1;
  } else {
    xl=i1; xh=i0;
  }
  /*initial guess for NR:*/
  x=0.5*(i0+i1); dxold=i1-i0; dx=dxold;
  p=lam*x*sin(x)+H*cos(x); dp=(lam-H)*sin(x)+lam*x*cos(x);
  /*find a root using N-R:*/
  for (its=1;its<=MAXIT;its++) {
    if ( ( ((x-xh)*dp-p)*((x-xl)*dp-p) >= 0.0 )
	|| ( fabs(2.0*p) > fabs(dxold*dp) ) ) { 
	/*bisect if NR jumps out of range, or if not decreasing fast enough.*/
	dxold=dx; dx=0.5*(xh-xl);
	x=xl+dx;
	if (x == xl) {w[0] = x; return;}/*change in root is negligible.*/
    } else {
	dxold=dx; dx=p/dp;
	x0=x; x=x0-dx;
	if (x == x0) {w[0]=x; return;}/*change in root is negligible.*/
    }
    if (fabs(dx) < EPS) {w[0]=x; return;}
    p=lam*x*sin(x)+H*cos(x); dp=(lam-H)*sin(x)+lam*x*cos(x);
    if (p < 0.0 ) {/*maintain the bracket on the root*/
	xl=x;
    } else {
	xh=x;
    }
  }
  w[0]=0.0; return; /*has done MAXIT runs but hasn't converged yet*/
}

/*define inputs and outputs*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  double *w;
  lam = *mxGetPr(prhs[0]);
  H = *mxGetPr(prhs[1]);
  i0 = *mxGetPr(prhs[2]);
  i1 = *mxGetPr(prhs[3]);
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  w = mxGetPr(plhs[0]);

  imag_root_wtr(w,lam,H,i0,i1);

  return;
}
