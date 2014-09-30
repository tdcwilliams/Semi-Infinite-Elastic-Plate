% this module provides the subroutine GEN_gauleg, which computes
% N-vectors x & w of abscissae and weights
% to use for Gauss-Legendre integration
% given N and the endpoints x1 & x2.
% CALL [x,w]=GEN_gauleg(x1,x2,N)function

function [x,w]=GEN_numint_chebyshev(N)
%% \int_{-1}^1f(x)/\sqrt{1-x^2}dx~w'*f(x)
%% source of formula:
%% http://www.efunda.com/math/num_integration/num_int_gauss.cfm

disp('please use OP_numint_chebyshev.m vs GEN_*')

if N==0
  x=0;
  w=pi;
  return;
end

nvec=(1:N)';
x=cos( pi/2/N*(2*nvec-1) );
w=0*x+pi/N;

if 0
  cos(N*acos(x))
end

return