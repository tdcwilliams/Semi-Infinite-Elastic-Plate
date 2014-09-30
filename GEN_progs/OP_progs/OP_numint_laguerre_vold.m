function [x,w]=OP_numint_laguerre_vold(alpha,n)
%% CALL: [x,w]=OP_numint_laguerre(alpha,n)
%% \int_0^\infty.f(\xi).xi^alpha*e^{-\xi}.d\xi ~ w'*f(x).

ip=-prod(1+alpha./(1:n-1))*gamma(alpha+1)/n;
[x,w0]=numint_gausslaguerre0(alpha,n);
w=ip*w0;

function [x,w]=numint_gausslaguerre0(alpha,n);
