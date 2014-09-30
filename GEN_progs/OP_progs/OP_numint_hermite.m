function [x, w] = OP_numint_hermite(n)
%% CALL: [xx, ww] = OP_numint_hermite(alpha,n)
%%  Gaussian quadrature routine that approx's the
%%   integral \int_{-\infty}^\infty.w(t).f(t).dt ~ ww'*f(xx),
%%   where the weight function is w(t)=e^{-t^2};
%%  n>=2 is no of points;
%% ########################################################
%% code posted by Geert V.D. on 28/6/08 on
%% http://www.mathworks.com/matlabcentral/fileexchange/8067
%% ########################################################

i = 1:n-1;
a = sqrt(i/2);
CM = diag(a,1) + diag(a,-1);

[V L] = eig(CM);
[x ind] = sort(diag(L));
V = V(:,ind)';
w = sqrt(pi) * V(:,1).^2;