function [xbig,xsmall]=GEN_splitorders(x);

DO_TEST=0;
if nargin==0
  x=pi/1e13
  DO_TEST=1;
end

expon=floor(log10(x));
xbig=floor(x.*10.^(7-expon))*10.^round(expon-7)
xsmall=x-xbig

if DO_TEST
  x1=x+1e-22
  [x1b,x1s]=GEN_splitorders(x1);
  x1-x
  {x1b-xbig,x1s-xsmall}
%  get(x1s)
end