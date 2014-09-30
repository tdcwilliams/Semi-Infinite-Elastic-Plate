function [x,w]=OP_numint_chebyshev(N)
%% CALL: [x,w]=OP_numint_chebyshev(N)
%% x = [roots of T_N(X)=0], w are weights so that
%% \int_{-1}^1f(x)/\sqrt{1-x^2}dx~w'*f(x)=\sum_{n=1}^N[w_n.f(x_n)];
%% source of formula:
%% http://www.efunda.com/math/num_integration/num_int_gauss.cfm

if N==0
  x   = 0;
  w   = pi;
  return;
end

nvec  = (N:-1:1)';
x     = cos( pi/2/N*(2*nvec-1) );
w     = 0*x+pi/N;

if 0
  cos(N*acos(x))
end

return
