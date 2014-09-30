function [A,f]=Test_fn2(y,varargin)

z=0*y;
f=[cos(y)-sin(y).^2, cos(y)-cos(y).*sin(y)];
A=[sin(y),z,z,cos(y)];

%% general solution is [c1*e^{-cos(x)};c2*e^{sin(x)}]