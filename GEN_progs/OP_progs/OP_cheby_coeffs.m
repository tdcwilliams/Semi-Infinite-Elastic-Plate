function [x,ipT,hnT,t]=OP_cheby_coeffs(a,b,Npolys,Nint)
%%
if nargin==3
  Nint=Npolys+1;
end

[t,w]=OP_numint_chebyshev(Nint);
[ipT,hnT]=OP_inprod_chebyshev(t,w,Npolys);
x=a+(b-a)/2*(1+t);