function x=GEN_findroots_muller(FXN,x3,varargin);

DO_TEST=0;
if nargin==0
  DO_TEST=1;
  FXN=@tester;
  x3=[[.32 .21+.001i,.6];[.32 .21+.001i,.6] ];
  prams=[];
end
prams=[];
%%
xk=x3(:,1);
xkm1=x3(:,2);
xkm2=x3(:,3);
%%
TOL=1e-13;
MAXITS=100;
critter=1;
its=0;
%%
while critter
  fk=feval(FXN,xk,varargin{:});
  fkm1=feval(FXN,xkm1,varargin{:});
  fkm2=feval(FXN,xkm2,varargin{:});
  %%
  dd01=(fkm1-fk)./(xkm1-xk);
  dd02=(fkm2-fk)./(xkm2-xk);
  dd12=(fkm2-fkm1)./(xkm2-xkm1);
  dd012=(dd12-dd01)./(xkm2-xk);
  %%
  ww=dd01+dd02-dd12;
  denom=ww+sqrt(ww.^2-4*fk.*dd012);
  denom_m=ww-sqrt(ww.^2-4*fk.*dd012);
  %%
  Jp=find(abs(denom_m)>abs(denom));
  denom(Jp)=denom_m(Jp);
  dx=-2*fk./denom;
  x=xk+dx;
  %%
  ERR=abs(dx);
  Jnot=find( abs(dx)>TOL | abs(dx./x)>TOL );
  critter=(~isempty(Jnot) & its<MAXITS);
  its=its+1;
  xkm2=xkm1;
  xkm1=xk;
  xk=x;
end

if its==MAXITS
  disp('warning (GEN_findroot_muller.m): root not converged');
%    [guess(unfinished_roots),x(unfinished_roots)]
  xnot=x(Jnot);
  x(Jnot)=NaN;
%  [guess,x]
end

function y=tester(zz,prams)
y=(zz-.3).^2;
