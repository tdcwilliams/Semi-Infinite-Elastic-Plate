function [xvec,exitflag]=...
  GEN_findroots_bisection(fxn,interval,varargin)

tol=1e-12;
%%
%z=varargin{:}%%extra arguments for function to be zeroed.
i0=interval(:,1);
i1=interval(:,2);
xvec=0*i1;
exitflag=xvec;

fl=feval(fxn,i0,varargin{:});
fh=feval(fxn,i1,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check there is a root in interval
critter=(fl > 0 & fh > 0) | (fl < 0 & fh < 0 );
J_NOTOK=find(critter);
J_notcvg=find( (fl.*fh <= 0) );

if ~isempty(J_NOTOK)
  xvec(J_NOTOK)=NaN;
  exitflag(J_NOTOK)=0;
  disp('warning (GEN_findroots_bisection.m):');
  disp('no root in intervals');
  disp([i0(J_NOTOK),i1(J_NOTOK)]);
  disp([fl(J_NOTOK),fh(J_NOTOK)]);pause
  %%
  fl(J_NOTOK)=[];
  fh(J_NOTOK)=[];
  i0(J_NOTOK)=[];
  i1(J_NOTOK)=[];
end

if isempty(i0)
  return;
elseif 0
  'there is a root (CHG OF SIGN)',pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if either endpoint is "close enough", just use that for root
J_i0zero=find(fl==0);
xvec(J_notcvg(J_i0zero))=i0(J_notcvg(J_i0zero));
exitflag(J_notcvg(J_i0zero))=1;
%%
J_i1zero=find(fh==0);
xvec(J_notcvg(J_i1zero))=i1(J_notcvg(J_i0zero));
exitflag(J_notcvg(J_i1zero))=2;
%%
fl([J_i0zero,J_i1zero])=[];
fh([J_i0zero,J_i1zero])=[];
J_notcvg([J_i0zero,J_i1zero])=[];
i0([J_i0zero,J_i1zero])=[];
i1([J_i0zero,J_i1zero])=[];

if isempty(i0)
  return;
elseif 0
  'there is a root (ENDS ARE NOT ZEROS)',pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define search direction to be in dirn of increasing p
xl=i0;
xh=i1;
xl(find(fl>0))=i1(find(fl>0));
xh(find(fl>0))=i0(find(fl>0));
x=0.5*(xl+xh);
dx=(xh-xl)/2;
%%
critter=1;%{x,dx,xl,xh}
Nits=0;
while critter
  Nits=Nits+1;
%    if Nits==36;{x,dx,xl,xh},end;
  p=feval(fxn,x,varargin{:});

  %% if p(x)==0 then have found the root;
  J_gszero=find(p==0);
  xvec(J_notcvg(J_gszero))=x(J_gszero);
  exitflag(J_notcvg(J_gszero))=3;
  %%
  xl(J_gszero)=[];
  xh(J_gszero)=[];
  x(J_gszero)=[];
  p(J_gszero)=[];
  J_notcvg(J_gszero)=[];
  dx(J_gszero)=[];

  if isempty(x)
    critter=0;
    break;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Jpos=find(p>0);
  xh(Jpos)=x(Jpos);
  %%
  Jneg=find(p<0);
  xl(Jneg)=x(Jneg);
  %%
  dx=(xh-xl)/2;%,pause
  x=xl+dx;
  %%
  Jcvg=find(abs(dx)<=tol);
  xvec(J_notcvg(Jcvg))=x(Jcvg);
  %%
  xl(Jcvg)=[];
  xh(Jcvg)=[];
  x(Jcvg)=[];
  dx(Jcvg)=[];
  %%
  J_notcvg(Jcvg)=[];
  critter=~isempty(x);
end
%  while abs(dx)>tol
%    p=feval(fxn,x,varargin{:});
%    if p==0
%      exitflag=3;
%      return;
%    elseif p>0
%      xh=x;
%    else
%      xl=x;
%    end
%    dx=(xh-xl)/2;
%    x=xl+dx;
%  end
%  exitflag=4;

%  xvec,pause