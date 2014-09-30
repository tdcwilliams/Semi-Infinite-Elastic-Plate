function [x, w] = OP_numint_laguerre(alpha,n)
%% CALL: [xx, ww] = OP_numint_laguerre(alpha,n)
%%  Gaussian quadrature routine that approx's the
%%   integral \int_0^\infty.w(t).f(t).dt ~ ww'*f(xx),
%%   where the weight function is w(t)=t^alpha*e^{-t};
%%  n>=2 is no of points;
%% ########################################################
%% code posted by Geert V.D. on 28/6/08 on
%% http://www.mathworks.com/matlabcentral/fileexchange/8067
%% ########################################################

n=max(2,n);%% n>=2

i = 1:n;
a = (2*i-1) + alpha;
b = sqrt( i(1:n-1) .* ((1:n-1) + alpha) );
CM = diag(a) + diag(b,1) + diag(b,-1);

[V L] = eig(CM);
[x ind] = sort(diag(L));
V = V(:,ind)';
w = gamma(alpha+1) .* V(:,1).^2;