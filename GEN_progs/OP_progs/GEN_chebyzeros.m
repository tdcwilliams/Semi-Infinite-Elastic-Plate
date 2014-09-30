function y=GEN_chebyzeros(a,b,N);
%% CALL: y=GEN_chebyzeros(a,b,N);
%% when using N chebyshev polynomials to interpolate a function f
%% over an interval [a,b],
%% this function returns the points inside [a,b] at which f must be evaluated.

cheby_zeros=cos( pi/2/N*(2*(1:N)'-1) );
y=a+(b-a)*(cheby_zeros+1)/2;
