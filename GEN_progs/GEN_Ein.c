/* To Compile:
 * 'mex Gaulag0.c' (on a machine with gsl installed) */
#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>

#define EPS 3.0e-14
#define MAXIT 5000

int M,N;

void Ein(double *,double *,double *,double *,double *,double *);

void Ein(double *yyr, double *yyi, double *zzr, double *zzi, double *testr, double *testi) {

int i, its, j, r;
double yr ,yi, zr, zi, dyr, dyi, dyr1;

 for (i=1;i<=M;i++) {
   for (j=1;j<=N;j++) {
     r=i-1+M*(j-1);
     zr=zzr[r]; zi=zzi[r];
     dyr=zr; dyi=zi;
     yr=dyr; yi=dyi;
     testr[0]=dyr; testi[0]=dyi;
     for (its=2;its<=MAXIT;its++) {
	/*dy=-(its-1)*1.0/(its*its)*z*dy;*/
	dyr1=-(its-1)*1.0/(its*its)*(zr*dyr-zi*dyi);
	dyi=-(its-1)*1.0/(its*its)*(zr*dyi+zi*dyr); dyr=dyr1;
	yr=yr+dyr; yi=yi+dyi;
	testr[its-1]=dyr; testi[its-1]=dyi;
	/*testr[its-1+MAXIT]=zr;
	testr[its-1+2*MAXIT]=zr;
	testr[its-1+3*MAXIT]=;
	testr[its-1+4*MAXIT]=zr;
	testr[its-1+5*MAXIT]=zr;*/
	if(sqrt(dyr*dyr+dyi*dyi)<EPS) break;
     }
     yyr[r]=yr; yyi[r]=yi;
   }
 }
  return;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{
  double *yyr, *yyi, *zzr, *zzi, *testr, *testi;
  zzr=mxGetPr(prhs[0]);
  zzi=mxGetPi(prhs[0]);
  M=mxGetM(prhs[0]);
  N=mxGetN(prhs[0]);

  plhs[0]=mxCreateDoubleMatrix(M,N,mxCOMPLEX);
  yyr = mxGetPr(plhs[0]);
  yyi = mxGetPi(plhs[0]);
  plhs[1]=mxCreateDoubleMatrix(MAXIT,1,mxCOMPLEX);
  testr = mxGetPr(plhs[1]);
  testi = mxGetPi(plhs[1]);

  Ein(yyr, yyi, zzr, zzi, testr, testi);

  return;
}
