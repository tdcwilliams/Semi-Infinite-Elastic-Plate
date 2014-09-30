function [xvec,exitflag]=...
           GEN_findroot_NRsafe2_num(fxn,interval,guess,target,varargin)
%% same as GEN_findroot_NRsafe.m but finds multiple roots simultaneously;
%% also if fxn=@f, this program solves f(xvec)==target
%%  (target a vector of constants),
%%  instead of just f(xvec)==0.

%% CONVERGENCE PARAMETERS:
EPS=1e-12;
MAXIT=300;
Nroots=length(guess);
xvec=0*guess;
exitflag=xvec;

i0=interval(1);
i1=interval(2);
%%
fl=feval(fxn,i0,varargin{:})-target;
fh=feval(fxn,i1,varargin{:})-target;

%% check there is a root in interval
J_NOTOK=find( (fl.*fh > 0) );
J_OK=find( (fl.*fh <= 0) );
guess(J_NOTOK)=[];
fl(J_NOTOK)=[];
fh(J_NOTOK)=[];
target(J_NOTOK)=[];
%%
xvec(J_NOTOK)=NaN;
%% exitflag=0 in this case, so no need to change;

%% if either endpoint is "close enough", just use that for root
J_i0zero=find(fl==0);
xvec(J_OK(J_i0zero))=i0;
exitflag(J_OK(J_i0zero))=1;

J_i1zero=find(fh==0);
xvec(J_OK(J_i1zero))=i1;
exitflag(J_OK(J_i1zero))=2;
%%
guess([J_i0zero,J_i1zero])=[];
fl([J_i0zero,J_i1zero])=[];
fh([J_i0zero,J_i1zero])=[];
J_OK([J_i0zero,J_i1zero])=[];
target([J_i0zero,J_i1zero])=[];
if isempty(guess)
  return;
end

%% define search direction to be in dirn of increasing p
xl=i0+0*guess;
xh=i1+0*guess;
xl(find(fl>0))=i1;
xh(find(fl>0))=i0;

%% initial guess for NR:
x=guess;
p=feval(fxn,x,varargin{:});
p2=feval(fxn,x+EPS/2,varargin{:});
p1=feval(fxn,x-EPS/2,varargin{:});
dp=(p2-p1)/EPS;
p=p-target;
%%
J_gszero=find(p==0);
xvec(J_OK(J_gszero))=x(J_gszero);
exitflag(J_OK(J_gszero))=3;
%%
xl(J_gszero)=[];
xh(J_gszero)=[];
x(J_gszero)=[];
p(J_gszero)=[];
dp(J_gszero)=[];
J_OK(J_gszero)=[];
target(J_gszero)=[];
if isempty(x)
  return;
end
%%
dxold=i1-i0+0*x;
dx=dxold;

%% find a root using N-R:*/
for its=1:MAXIT
  crit1=( ((x-xh).*dp-p).*((x-xl).*dp-p) >= 0 );
  crit2=( abs(2*p) > abs(dxold.*dp) );
  %%
  J_bisect=find( crit1|crit2 );
  J_NR=find( ~crit1 & ~crit2 );
  %%
  dxold=dx;
  dx(J_bisect)=0.5*(xh(J_bisect)-xl(J_bisect));
  x(J_bisect)=xl(J_bisect)+dx(J_bisect);
  %%
  dx(J_NR)=p(J_NR)./dp(J_NR);
  x0=x(J_NR);
  x(J_NR)=x0-dx(J_NR);

  %% TEST FOR CONVERGENCE:
  J_cvg=find(abs(dx) < EPS);
  xvec(J_OK(J_cvg))=x(J_cvg);
  exitflag(J_OK(J_cvg))=3;
  %%
  xl(J_cvg)=[];
  xh(J_cvg)=[];
  x(J_cvg)=[];
  J_OK(J_cvg)=[];
  target(J_cvg)=[];
  dx(J_cvg)=[];
  dxold(J_cvg)=[];
  if isempty(x)
    return;
  end

  p=feval(fxn,x,varargin{:});
  p2=feval(fxn,x+EPS/2,varargin{:});
  p1=feval(fxn,x-EPS/2,varargin{:});
  dp=(p2-p1)/EPS;
  p=p-target;
  %%
  J_gszero=find(p==0);
  xvec(J_OK(J_gszero))=x(J_gszero);
  exitflag(J_OK(J_gszero))=3;
  %%
  xl(J_gszero)=[];
  xh(J_gszero)=[];
  x(J_gszero)=[];
  p(J_gszero)=[];
  dp(J_gszero)=[];
  J_OK(J_gszero)=[];
  target(J_gszero)=[];
  dx(J_gszero)=[];
  dxold(J_gszero)=[];
  if isempty(x)
    return;
  end

  %% maintain the bracket on the root
  J0=find(p<0);
  xl(J0)=x(J0);
  J0=find(p>0);
  xh(J0)=x(J0);
end

%% if program reaches here it
%% has done MAXIT runs but hasn't converged yet
xvec(J_OK)=NaN;
exitflag(J_OK)=4;
disp('warning (GEN_findroot_NR2.m): some roots not converged');