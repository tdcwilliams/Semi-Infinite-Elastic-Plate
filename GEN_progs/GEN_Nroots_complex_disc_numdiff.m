function Nroots=...
GEN_Nroots_complex_disc_numdiff(FXN,centre,radius,varargin);
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
tol=1e-3;
maxN=500;
%%
dN=10;
N=0;
err=2*tol;
Nroots=2e7;
%%
M=4;
Nroots_old=Nroots+(0:M)';
%%
for j=0:M
  N=N+(j+1)*dN;
  jj=N+1:2*N+1;
  [t,w]=GEN_numint_exp(N+1);
%    ip=GEN_inprod_exp(t,w);
  %%
  kk=centre+radius*exp(i*pi*t);
  ff=feval(FXN,kk,varargin{:});
  %%
  dk=0*kk;
  dk(2:end-1)=kk(3:end)-kk(1:end-2);
  dk(1)=kk(2)-kk(end);
  dk(end)=kk(1)-kk(end-1);
  %%
  df=0*ff;
  df(2:end-1)=ff(3:end)-ff(1:end-2);
  df(1)=ff(2)-ff(end);
  df(end)=ff(1)-ff(end-1);
  df=df./dk;
  %%
  integral=radius/2*sum(w.*exp(i*pi*t).*df./ff);
  Nroots_old(M+1-j)=round(real( integral ));
end
critter=1;

while critter & (N<maxN)
  N=N+dN;
  jj=N+1:2*N+1;
  [t,w]=GEN_numint_exp(N+1);
%    ip=GEN_inprod_exp(t,w);
  %%
  kk=centre+radius*exp(i*pi*t);
  ff=feval(FXN,kk,varargin{:});
  %%
  dk=0*kk;
  dk(2:end-1)=kk(3:end)-kk(1:end-2);
  dk(1)=kk(2)-kk(end);
  dk(end)=kk(1)-kk(end-1);
  %%
  df=0*ff;
  df(2:end-1)=ff(3:end)-ff(1:end-2);
  df(1)=ff(2)-ff(end);
  df(end)=ff(1)-ff(end-1);
  df=df./dk;
  %%
  integral=radius/2*sum(w.*exp(i*pi*t).*df./ff);
  Nroots=round(real( integral ));
  Nroots_old=[Nroots;Nroots_old(1:M)];
  critter=~( sum(Nroots_old==Nroots)==(M+1) );
end

if N==maxN
  Nroots=NaN;
end