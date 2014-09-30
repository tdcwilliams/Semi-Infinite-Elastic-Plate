function x=...
GEN_findroot_complex_disc(FXN,centre,radius,Niter,varargin);
%% 1. approx's fxn on a disc using Cauchy integral formula
%%   to generate a polynomial;
%% 2. uses roots to find zeros & discards those outside
%%   the disc.
%% NB it is able to find triple roots
%%  but not quadruple roots;
%% NB can't resolve roots less than 1e-12 apart

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
tol=1e-12;
if radius<=tol
  Nroots=GEN_Nroots_complex_disc_numdiff(FXN,...
    centre,radius,varargin{:});
  x(1:Nroots)=centre;
  return;
end

septol=1e-3;
%% this helps to determine if roots are distinct
%%  or if there are multiple roots;
%%   septol=1e-6 is able to find double roots;
%%   septol=1e-3 is able to find triple roots;
maxNiter=10;%Niter
maxN=200;
%%
dN=10;
N=0;
err=2*tol;
while err>tol & N<maxN
  N=N+dN;
  jj=N+1:2*N+1;
  [t,w]=GEN_numint_exp(N+1);
  ip=GEN_inprod_exp(t,w);
  %%
  kk=centre+radius*exp(i*pi*t);
  ff=feval(FXN,kk,varargin{:});
  fn_tilde=ip*ff;
  %%
%    plot(t,real(ff)), hold on;
%    plot(t,real(GEN_interp_exp(t,fn_tilde)),'--r');
%    pause
  %%
  fn_tilde=fn_tilde(jj);
  err=abs(max(fn_tilde(end-1:end))/max(fn_tilde));
end

N=max(find(fn_tilde~=0))-1;
fn_tilde=fn_tilde(1:N+1);

if 0
  CM=[ [zeros(1,N-1);eye(N-1)],...
        -fn_tilde(1:N)/fn_tilde(N+1) ];
  xtilde=eig(CM);
else
  xtilde=roots(flipud(fn_tilde));
end
xtilde(find(abs(xtilde)>1))=[];
Nrts=length(xtilde);
x=centre+radius*xtilde;%,pause

%% test if roots are resolved correctly
[Xtilde,Xtilde0]=meshgrid(xtilde,xtilde);
R=abs(Xtilde-Xtilde0);
%  R<septol & R>0
[j1,r1]=find(R<septol);
JP=find(r1>j1);
j1=j1(JP);
r1=r1(JP);
%%
if 0%N==maxN
  disp('warning (GEN_findroot_complex_disc.m): Fourier series not converged');x,radius,pause
%    x=GEN_findroot_NRnum(FXN,x,varargin{:});
  for j=1:length(x)
    x(j)=GEN_findroot_complex_disc(FXN,x(j),...
           radius/10,Niter+1,varargin{:});pause
  end
end


if isempty(JP)
  return;
end
%%
S=1;
s=1;
XNR{S}=[x(j1(s)),x(r1(s))];
JR{S}=[j1(s),r1(s)];

for s=2:length(r1)
  jr=[j1(1:s-1) r1(1:s-1)];
%    jc=[j1(s) r1(s)]
  crit1=isempty(find(j1(s)==jr));
  crit2=isempty(find(r1(s)==jr));
  if (crit1&crit2)%% both roots completely new
    S=S+1;
    XNR{S}=[x(j1(s)),x(r1(s))];
    JR{S}=[j1(s),r1(s)];
  elseif crit1%% only 'j1' root is new
    XNR{S}=[XNR{S},x(j1(s))];
    JR{S}=[JR{S},j1(s)];
  elseif crit2%% only 'r1' root is new
    XNR{S}=[XNR{S},x(r1(s))];
    JR{S}=[JR{S},r1(s)];
  end
end

%  for j=1:S
%    [j S]
%    XNR{j}
%    JR{j}
%  end,pause

%% NOW TEST TO SEE HOW CLOSE TOGETHER ROOTS ARE - THIS IS TO
%%  RESOLVE MULTIPLE ROOTS:
for j=1:S
  xnr=XNR{j};
  xav=mean(xnr);
  rad=abs(xnr(1)-xav);%,[2*rad,tol]
  critter=(2*rad<tol);%%roots close enough together
%    Niter
  if critter;
    x(JR{j})=xav;
  else

tst=(abs(xnr-xav)<10*rad)
    y=GEN_findroot_complex_disc(FXN,xav,2*rad,...
        Niter+1,varargin{:});
    x(JR{j})=y;
  end
end

%  disp('START')
%  for j=1:Nrts
%    x(j)
%  end