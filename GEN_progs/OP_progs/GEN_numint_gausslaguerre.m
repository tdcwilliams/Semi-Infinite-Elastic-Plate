function [x,w]=GEN_numint_gausslaguerre(alpha,n)
%% CALL: [x,w]=GEN_numint_gausslaguerre(alpha,n)
%% \int_0^\infty.f(\xi).xi^alpha*e^{-\xi}.d\xi ~ w'*f(x).

disp('please use OP_numint_laguerre.m')

ip=-prod(1+alpha./(1:n-1))*gamma(alpha+1)/n;
[x,w0]=Gaulag0(alpha,n); w=ip*w0;
