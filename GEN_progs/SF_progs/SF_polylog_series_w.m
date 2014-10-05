%% SF_polylog_series_w.m
%% Author: Timothy Williams
%% Date:   20140913, 22:13:20 CEST

function [y,k]=SF_polylog_series_w(z,n,OPT);
%% algorithms from Wood (1992)
%% OPT==1 (for z~1):  series (9.5) for  n positive integers, else (9.3)
%% OPT==0 (for z~-1): series (9.2)
%% TODO: check branch cut OK for z~1

tol   = 1e-12;

if nargin==0
   if 0
      z     = exp(1i*pi/4)*[1;1];
      n     = 2.5;
      OPT   = 1;
   elseif 1
      Hr    = 1.011010578882925;
      %Hr    = 1.000000578882925;
      z     = exp(2i*pi/Hr);%%slow if near 1;
      %n     = 7/3;
      n     = 3
      OPT   = 1
   else
      Hr    = 1.011010578882925;
      %Hr    = 1.000000578882925;
      z     = -exp(2i*pi/Hr);%%slow if near 1;
      n     = 7/3;
      OPT   = 0;
   end
end

k     = 0;
t1    = 1;
if OPT==1%%good for near 1
   x        = log(z);%%calculate Li_n(e^x)
   POS_INT  = ((n==round(n))&(n>0));
   if POS_INT==0
      y  = zeta(n)+gamma(1-n)*(-x).^(n-1);
   else
      t0 = 1;
      H  = 0;
      for k=1:n-1
         t0 = t0.*x/k;%%x^k/k!
         H  = H+1/k;
      end
      y  = zeta(n)+(H-log(-x)).*t0;
   end
else%%good for near -1
   x  = log(-z);%%calculate Li_n(e^x)
   y  = -fn_eta(n);
end
%%
critter  = 1;

while critter
   k  = k+1;
   t1 = (x/k).*t1;%%x^k/k!
   if OPT==1
      if POS_INT==0
         dy = zeta(n-k)*t1;
         y  = y+dy;
      elseif k~=(n-1)
         dy = zeta(n-k)*t1;
         y  = y+dy;
      end
   else
      dy = -fn_eta(n-k)*t1;
      y  = y+dy;
   end
   %%
   critter  = (max(abs(dy))>tol);
   % {k,max(abs(dy)),y(1)}
end

function y=fn_eta(p)

y  = (1-2^(1-p))*zeta(p);
