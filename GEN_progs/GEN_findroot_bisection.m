function [x,exitflag]=...
GEN_findroot_bisection(fxn,interval,varargin)

tol   = 1e-12;
%%
%z=varargin{:}%%extra arguments for function to be zeroed.
i0 = interval(1);
i1 = interval(2);

fl = feval(fxn,i0,varargin{:});
fh = feval(fxn,i1,varargin{:});

critter  = ((fl > 0 & fh > 0) | (fl < 0 & fh < 0 ));
tries    = 0;

%% check there is a root in interval
SHOW  = 0;
while critter
   tries = tries+1;
   if SHOW
      disp('trying to fix interval...');
      tries,interval,fvals = [fl,fh]
   end
   interval = fix_int(interval,[fl,fh]);
   %%
   i0 = interval(1);
   i1 = interval(2);
   fl = feval(fxn,i0,varargin{:});
   fh = feval(fxn,i1,varargin{:});

   critter  = ((fl > 0 & fh > 0) | (fl < 0 & fh < 0 ));

   %% check there is a root in interval
   if ~critter
      if SHOW
         disp('interval fixed!')
      end
   elseif tries==15
      x         = NaN;
      exitflag  = 0;

      disp('warning (GEN_findroot_bisection.m):');
      disp('no root in interval');
      interval,fvals = [fl,fh]
   
      return;
   end
   %pause
end

%% if either endpoint is "close enough", just use that for root
if (fl == 0)
  x         = i0;
  exitflag  = 1;
  return;
end
if (fh == 0)
  x         = i1;
  exitflag  = 2;
  return;
end

%% define search direction to be in dirn of increasing p
if (fl < 0.0)
  xl  = i0;
  xh  = i1;
else
  xl  = i1;
  xh  = i0;
end

x  = 0.5*(xl+xh);
dx = (xh-xl)/2;
%  {x,dx,xl,xh}

Nits     = 0;
critter  = 1;

while critter
  Nits   = Nits+1;
%    if Nits==36;{x,dx,xl,xh},end;
  p   = feval(fxn,x,varargin{:});
  if p==0
    exitflag   = 3;
    return;
  elseif p>0
    xh   = x;
  else
    xl   = x;
  end
  dx        = (xh-xl)/2;
  ERR       = abs(dx/x);
  critter   = ( abs(dx)>tol & ERR>tol );
  x         = xl+dx;%{x,dx,xl,xh}
end
exitflag = 4;

function int2   = fix_int(int,fvals)

int2  = int;
dx    = int(2)-int(1);
dy    = fvals(2)-fvals(1);
slope = dy/dx;
x0    = int(1)-fvals(1)/slope;

frac  = .1;
if x0<int(1)
   int2(1)  = x0-frac*dx;
   int2(2)  = int(1);
elseif x0>int(2)
   int2(2)  = x0+frac*dx;
   int2(1)  = int(2);
end
