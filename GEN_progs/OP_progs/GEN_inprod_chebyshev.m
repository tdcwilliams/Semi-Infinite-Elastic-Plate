function [An,hn]=GEN_inprod_chebyshev(xj,wj,NgP,fj)

disp('WARNING: GEN_inprod_chebyshev.m')
disp('please use OP_inprod_chebyshev.m')

if nargin==4
  [An,hn]=OP_inprod_chebyshev(xj,wj,NgP,fj);
else
  [An,hn]=OP_inprod_chebyshev(xj,wj,NgP);
end