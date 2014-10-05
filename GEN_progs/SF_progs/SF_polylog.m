function y=SF_polylog(x,n);
%% CALL: y=SF_polylog(x,n);
%%y=Li_n(x), n>0
%%
%% series:
%%  y=\sum_{k=1}^\infty {x^k/k^n = x^{k-1}/(k-1)^n*[((k-1)/k)^n*x]
%%   n=5 takes about       250 terms for x on the unit circle C;
%%   n=4 takes about     1 000 terms for x on the unit circle C;
%%   n=3 takes about    10 000 terms for x on the unit circle C;
%%   n=2 takes about 1 000 000 terms for x on the unit circle C;
%%
%% calculation with Bose integral is faster though;
%% y=x/gamma(n)*\int_0^\infty[t^{n-1}*exp(-t)/(1-x*exp(-t))]dt
%%
%% SF_polylog.m
%% Author: Timothy Williams
%% Date:   20140913, 22:13:20 CEST

tol   = 1e-12;

if nargin==0
   x  = [0;exp(1i*pi/4)*[1;1];1];
   n  = 1.5;
end

if n==1
   y  = -log(1-x);
   return;
end

J  = (1:length(x))';
y  = 0*x;

%%test for at 0:
j0    = find(x==0);
y(j0) = 0;
J(j0) = [];
x(j0) = [];
if isempty(x)
   return;
end

%%test for at 1:
j0       = find(x==1);
y(J(j0)) = zeta(n);
J(j0)    = [];
x(j0)    = [];
if isempty(x)
   return;
end

%%test for near 1:
tol1     = 1e-1;%%abs(log(z))<2*pi
j0       = find(abs(x-1)<tol1);
y(J(j0)) = SF_polylog_series_w(x(j0),n,1);
J(j0)    = [];
x(j0)    = [];
if isempty(x)
   return;
end

%%test for near -1:
tol1     = 1e-1;%%abs(log(-z))<pi
j0       = find(abs(x+1)<tol1);
y(J(j0)) = SF_polylog_series_w(x(j0),n,0);
J(j0)    = [];
x(j0)    = [];
if isempty(x)
   return;
end

if 0
   %%series expansion;
   jser     = find(abs(x<1)&(x~=0));
   y(jser)  = SF_polylog_series_z(x(jser),n);

   %%TODO: implement calculation for abs(x>0);
   j2    = find(abs(x)>1);
   y(j2) = NaN;
else
   %%Bose integral:
   jint     = find( (x~=0)&(x~=1)&(angle(x)~=0) );
   y(jint)  = SF_polylog_bose(x(jint),n);
end
