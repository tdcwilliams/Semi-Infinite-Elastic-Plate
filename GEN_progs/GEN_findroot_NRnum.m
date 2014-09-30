function [x,xnot,unfinished_roots]=...
  GEN_findroot_NRnum(fxn,guess,varargin)
%% finds roots with Newton-Rhapson, doing the differentiation numerically
%% CALL: xvec=GEN_findroots_NRnum(FXN,guess,v2,v3,...);
%% INPUTS:
%% FXN is a function handle, eg @f, where f=f(x,v2,v3,...);
%% guess is a guess for the root. It can be a column vector,
%%  in which case it finds multiple zeros.
%% v1,v2,v3... enter extra arguments of f after x3;
%% x is the root of 'fxn' nearest to 'guess'
%% NB 'guess' can be a vector
%%
%% OUTPUTS:
%% x = converged roots;
%% xnot = unconverged 'roots';
%% unfinished_roots=[] => all have converged
%%  unfinished_roots~=[] => some haven't converged -
%%  xnot correspond to guess(unfinished_roots);

if ~iscell(guess)
  TOL=1e-15;
else
  TOL=guess{2};
  guess=guess{1};
end
EPS=max(TOL/10,1e-13);
DX=EPS*exp(i*angle(guess));
%% increment for numerical differentiation

MAXITS=300;%% need this to stop infinite loops
unfinished_roots=1;
j=0;
x0=guess;

while ~isempty(unfinished_roots) & j<MAXITS
  p=feval(fxn,x0,varargin{:});
  %%
  p3=feval(fxn,x0+DX/2,varargin{:});
  p2=feval(fxn,x0-DX/2,varargin{:});
%    p4=feval(fxn,x0+DX,varargin{:});
%    p1=feval(fxn,x0-DX,varargin{:});
%    dp=(p4-8*p3+8*p2-p1)./(-12*DX);
  dp=(p3-p2)./DX;
  %%
  dx=-p./dp;
  x=x0+dx;
  %%
  j=j+1;
  x0=x;
  unfinished_roots=find(abs(dx)>TOL);
end

%  [j MAXITS]
if j==MAXITS
  disp('warning (GEN_findroot_NRnum.m): root not converged');
%    [guess(unfinished_roots),x(unfinished_roots)]
  xnot=x(unfinished_roots);
  x(unfinished_roots)=NaN;
%  [guess,x]
end
