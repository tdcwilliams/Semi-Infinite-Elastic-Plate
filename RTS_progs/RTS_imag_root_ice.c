/* To Compile:
 * 'mex imag_root_ice.c' (on a machine with gsl installed) */
/*CALL: x=imag_root_ice(del,H,i0,i1)
	->finds the imaginary root of the dispersion relation for ice,
	  f=del*w*sin(w)+H^5*cos(w), inside the interval [i0,i1]*/
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979
#define EPS 1.0e-8
#define MAXIT 100

/*declare input variables*/
double del, H, i0, i1;

/*construct main function*/
void imag_root_ice(double *, double, double, double, double, double *);

void imag_root_ice(double *w, double del, double H, double i0, double i1, double *exitflag) {

int its, j;
double p, dp, x, x0;
double xl, fl, xh, fh;
double dxold, dx, Lam, x4, H4;

  H4=H*H*H*H;
  x4=i0*i0*i0*i0;
  Lam=x4+del*H4;
  fl=Lam*i0*sin(i0)+H*H4*cos(i0);
  x4=i1*i1*i1*i1;
  Lam=x4+del*H4;
  fh=Lam*i1*sin(i1)+H4*H*cos(i1);
if (1) {/*necessary? may save time if either endpoint is "close enough"*/
  /*check there is a root in interval*/
  if ( (fl > 0.0 && fh > 0.0) || ( fl < 0.0) && fh < 0.0 ) {
	w[0]=0.0; exitflag[0]=0; /*exitflag[1]=fl; exitflag[2]=fh;*/ return;
  } 
  /*if either endpoint is "close enough", just use that for root*/
  if (fl == 0.0) {w[0] = i0; exitflag[0]=1; return;}
  if (fh == 0.0) {w[0] = i1; exitflag[0]=2; return;}
}
  /*define search direction to be in dirn of increasing p*/
  if (fl < 0.0) {
    xl=i0;
    xh=i1;
  } else {
    xl=i1;
    xh=i0;
  }
  /*initial guess for NR:*/
  x=0.5*(i0+i1);
  dxold=i1-i0;
  dx=dxold;
  x4=x*x*x*x;
  Lam=x4+del*H4;
  p=Lam*x*sin(x)+H4*H*cos(x);
  dp=(5*x4+H4*(del-H))*sin(x)+Lam*x*cos(x);
  /*find a root using N-R:*/
  for (its=1;its<=MAXIT;its++) {
    if ( ( ((x-xh)*dp-p)*((x-xl)*dp-p) >= 0.0 )
	|| ( fabs(2.0*p) > fabs(dxold*dp) ) ) { 
	/*bisect if NR jumps out of range, or if not decreasing fast enough.*/
	dxold=dx; dx=0.5*(xh-xl);
	x=xl+dx;
	if (x == xl) {w[0] = x; exitflag[0]=3; return;}/*change in root is negligible.*/
    } else {
	dxold=dx; dx=p/dp;
	x0=x; x=x0-dx;
	if (x == x0) {w[0]=x; exitflag[0]=3; return;}/*change in root is negligible.*/
    }
    if (fabs(dx) < EPS) {w[0]=x; exitflag[0]=3; return;}
    x4=x*x*x*x;
    Lam=x4+del*H4;
    p=Lam*x*sin(x)+H4*H*cos(x);
    dp=(5*x4+H4*(del-H))*sin(x)+Lam*x*cos(x);
    if (p < 0.0 ) {/*maintain the bracket on the root*/
	xl=x;
    } else {
	xh=x;
    }
  }
  w[0]=0.0;
  exitflag[0]=4;
  return; /*has done MAXIT runs but hasn't converged yet*/
}

/*define inputs and outputs*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  double *w;
  double *exitflag;/**/
  del = *mxGetPr(prhs[0]);
  H = *mxGetPr(prhs[1]);
  i0 = *mxGetPr(prhs[2]);
  i1 = *mxGetPr(prhs[3]);
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  w = mxGetPr(plhs[0]);
  plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
  exitflag = mxGetPr(plhs[1]);
/* 0 => no roots in i/val, 1 => w=i0, 2 => w=i1, 3 => w found normally,
	4 => MAXIT reached.*/

  imag_root_ice(w,del,H,i0,i1,exitflag);

  return;
}
