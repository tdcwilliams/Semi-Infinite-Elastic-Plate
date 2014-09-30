function [x,exitflag]=...
           GEN_findroot_NRsafe(fxn,interval,varargin)

%% CONVERGENCE PARAMETERS:
EPS=1e-12;
MAXIT=300;

i0=interval(1);
i1=interval(2);
%%
fl=feval(fxn,i0,varargin{:});
fh=feval(fxn,i1,varargin{:});
%%
%% check there is a root in interval
if (fl > 0 & fh > 0) | (fl < 0 & fh < 0 )
  x=0;
  exitflag=0;
  return;
end

%% if either endpoint is "close enough", just use that for root
if (fl == 0)
  x = i0;
  exitflag=1;
  return;
end
if (fh == 0)
  x = i1;
  exitflag=2;
  return;
end

%% define search direction to be in dirn of increasing p
if (fl < 0.0)
  xl=i0;
  xh=i1;
else
  xl=i1;
  xh=i0;
end

%% initial guess for NR:
x=0.5*(i0+i1);
[p,dp]=feval(fxn,x,varargin{:});

if p==0%%have found a root
  exitflag==3;
  return;
end
%%
dxold=i1-i0;
dx=dxold;

%% find a root using N-R:*/
for its=1:MAXIT
  crit1=( ((x-xh)*dp-p)*((x-xl)*dp-p) >= 0 );
  crit2=( abs(2*p) > abs(dxold*dp) );
  if crit1 | crit2
    %% bisect if NR jumps out of range,
    %% or if not decreasing fast enough.*/
    dxold=dx;
    dx=0.5*(xh-xl);
    x=xl+dx;
    if (x == xl)
      %% change in root is negligible.
      exitflag=3;
      return;
    end
  else
    %% use NR
    dxold=dx;
    dx=p/dp;
    x0=x;
    x=x0-dx;
  end
  if (abs(dx) < EPS)
    exitflag=3;
    return;
  end
  [p,dp]=feval(fxn,x,varargin{:});
  if p==0%%have found a root
    exitflag=3;
    return;
  end
  if (p < 0.0 ) %% maintain the bracket on the root
    xl=x;
  else
    xh=x;
  end
end
x=0;
exitflag=4;
disp('warning (GEN_findroot_NR.m): root not converged');
return;
%% if program reaches here it
%% has done MAXIT runs but hasn't converged yet