function xvec=GEN_findroots_muller(FXN,x3,varargin);
%% CALL: xvec=GEN_findroots_muller(FXN,x3,v2,v3,...);
%% INPUTS:
%% FXN is a function handle, eg @f, where f=f(x,v2,v3,...);
%% x3 is either [x1 x2 x3], 3 points near the root,
%%  or just x3, in which case x1 and x2 are created near x3.
%%  x3 can have multiple rows, in which case it finds multiple zeros.
%% v1,v2,v3... enter extra arguments of f after x3;

%  disp('yo, this is muller'),pause;
DO_TEST  = 0;
if nargin==0
  DO_TEST   = 1;
  FXN       = @tester;
  x3        = [[.32 .21+.001i,.6];[.32i .2i+.001,.6i] ];
  prams     = [];
end
prams = [];
%%
if size(x3,2)==1
   xk       = x3;
   xkm1     = x3*(1+1e-2);
   xkm2     = x3*(1-1e-3);
else
   xk       = x3(:,1);
   xkm1     = x3(:,2);
   xkm2     = x3(:,3);
end
%  guess = xk; [xk,xkm1,xkm2]
Jnot     = 1:length(xk);
xvec     = 0*xk;
%%
TOL      = 1e-13;
MAXITS   = 100;
critter  = 1;
its      = 0;
%%
while critter
  fk        = feval(FXN,xk,varargin{:});
  fkm1      = feval(FXN,xkm1,varargin{:});
  fkm2      = feval(FXN,xkm2,varargin{:});
  %%
  dd01      = (fkm1-fk)./(xkm1-xk);
  dd02      = (fkm2-fk)./(xkm2-xk);
  dd12      = (fkm2-fkm1)./(xkm2-xkm1);
  dd012     = (dd12-dd01)./(xkm2-xk);
  %%
  ww        = dd01+dd02-dd12;
  denom     = ww+sqrt(ww.^2-4*fk.*dd012);
  denom_m   = ww-sqrt(ww.^2-4*fk.*dd012);
  %%
  Jp        = find(abs(denom_m)>abs(denom));
  denom(Jp) = denom_m(Jp);
  dx        = -2*fk./denom;
%  guess
  x                  = xk+dx;
  %%
  ERR                = abs(dx);%[ERR,abs(dx./x)]
  Jcvg               = find( abs(dx)<TOL | abs(dx./x)<TOL );
  xvec(Jnot(Jcvg))   = x(Jcvg);%,pause
  %%
  Jnot(Jcvg)         = [];
  x(Jcvg)            = [];
  xk(Jcvg)           = [];
  xkm1(Jcvg)         = [];
  %%
  Jnan               = find(isnan(x));
  xvec(Jnot(Jnan))   = NaN;
  Jnot(Jnan)         = [];
  x(Jnan)            = [];
  xk(Jnan)           = [];
  xkm1(Jnan)         = [];
  %%
  its       = its+1;
  critter   = (~isempty(Jnot) & its<MAXITS);
  %%
  xkm2      = xkm1;
  xkm1      = xk;
  xk        = x;%[xk,xkm1,xkm2]
%  pause
%  xvec,pause
end

if its==MAXITS & ~isempty(Jnot)
  disp('warning (GEN_findroots_muller.m): root not converged');
%    [guess(unfinished_roots),x(unfinished_roots)]
  xnot         = x;
  xvec(Jnot)   = NaN;
%  [guess,x]
end

%% sometimes get NaN if guess is exactly on the root
%% so perturb the guess and try again
je = find(isnan(xvec));
if ~isempty(je)
   X3       = 1e-5+(1+1e-5)*x3(je,:);
   xvec(je) = GEN_findroots_muller(FXN,X3,varargin{:});
end


function y=tester(zz,prams)
y  = (zz-.3).^2.*(zz-.2i);
