%% SF_polylog_series_z.m
%% Author: Timothy Williams
%% Date:   20140913, 22:13:20 CEST

function y=SF_polylog_series_z(x,n);
%% series:
%%  y=\sum_{k=1}^\infty {x^k/k^n = x^{k-1}/(k-1)^n*[((k-1)/k)^n*x]
%%   n=5 takes about       250 terms for x on the unit circle C;
%%   n=4 takes about     1 000 terms for x on the unit circle C;
%%   n=3 takes about    10 000 terms for x on the unit circle C;
%%   n=2 takes about 1 000 000 terms for x on the unit circle C;

tol   = 1e-12;

if nargin==0
   x  = exp(1i*pi/4)*[1;1];
   n  = 3;
end

k        = 1;
y        = x;
dy       = x;
critter  = 1;

while critter
   k        = k+1;
   dy       = dy.*x*((k-1)/k)^n;
   y        = y+dy;
   critter  = (max(abs(dy))>tol);
   % {k,max(abs(dy)),y(1)}
end

%k,y
