%% SF_polylog_bose.m
%% Author: Timothy Williams
%% Date:   20140916
%%
%% Calculates y=polylog(x,n)=Li_n(x) with numerical (Gauss_Laguerre) integration
%% Uses the Bose integral representation:
%% y=x/gamma(n)*\int_0^\infty[t^{n-1}*exp(-t)/(1-x*exp(-t))]dt
%%
%% CALL: y=SF_polylog_bose(x,n,Nint);
%% y=Li_n(x);
%% (x>1 & x real) not permitted;
%% all other complex values of x are allowed;
%% Nint is the number of integration points;
%% if not specified, iterate to give sufficient convergence

function [y,Nint]=SF_polylog_bose(x,n,Nint);
tol   = 1e-12;

if nargin==0
   if 0
      x  = exp(1i*pi/4)*[1;1];
      n  = 1.5;
   else
      Hr = 1.011010578882925;
      x  = exp(2i*pi/Hr);%%slow if near 1;
      n  = 7/3;
   end
   %Nint  = 1000;%%tst converged
end

DO_CVG   = 1;
if exist('Nint','var')
   DO_CVG   = 0;
   Nint;
else
   %Nint  = 1000
   %Nint  = 200
   Nint  = 150;
end

fac   = x/gamma(n);

if 1%%straight Laguerre integration
   [t,w] = OP_numint_laguerre(n-1,Nint);
   I     = igrand_0(t,x)*w;
   y     = fac.*I;
end

critter  = 1;
dNint    = 50;
while critter&DO_CVG
   y0       = y;
   Nint     = Nint+dNint;
   y        = SF_polylog_bose(x,n,Nint);
   critter  = (max(abs(1-y0./y))>tol);
end

function y  = igrand_0(t,x)

[T,X] = meshgrid(t,x);
y     = 1./(1-X.*exp(-T));
