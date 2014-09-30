function Nroots=...
GEN_Nroots_complex_disc(FXN,centre,radius,varargin);
%% determines the number of roots inside a disc in the complex plane;
%% uses the formula N=1/(2\pi\rmi)\oint(f'/f)dz

if nargin==0%%use test arguments:
  FXN=@test_function_findroot;
  centre=0;
  radius=1.5*pi;
  Niter=0;
end

if nargin<4
  Niter=0;
end

%  Niter
tol=1e-6;
maxN=200;
%%
dN=10;
N=0;
err=2*tol;
Nroots0=1e7;
%%
while (err>tol || abs(imag(Nroots0))>1e-6) & (N<maxN)
  N=N+dN;
  jj=N+1:2*N+1;
  [t,w]=GEN_numint_exp(N+1);
%    ip=GEN_inprod_exp(t,w);
  %%
  kk=centre+radius*exp(i*pi*t);
  [ff,df]=feval(FXN,kk,varargin{:});
  %%
  Nroots=radius/2*sum(w.*exp(i*pi*t).*df./ff);
  err=min(abs(Nroots-Nroots0),abs(1-Nroots/Nroots0));
  Nroots0=Nroots;
end

if N<maxN
  Nroots=round(abs(Nroots));
else
  Nroots=NaN;
end