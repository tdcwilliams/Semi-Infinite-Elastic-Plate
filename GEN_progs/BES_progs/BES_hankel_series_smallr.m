function [Y,sing_coeff,err_report]=...
BES_hankel_series_smallr(r,A,gam,n,keepsing)
%% call: [Y,log_coeff]=BES_hankel_series(x,An,gam,m,keepsing)
%% calc's Y=Lw=(1/r.d/dr)^nW, where
%% W=sum_{n=-2}^length(gam).A_n.H_0(gam_n|x|)
%% uses the result that L.H_0(kr)=(-k/r)^nH_n(kr)
%% used to calc the thin plate 'in vacuo' Green's function
%% or the full floating thin plate Green's function at the surface.

TOL=1e-8;

sing_coeff=zeros(n+1,1);

jp=find(r>0);
rp=r(jp);
Y=zeros(size(r));
if isempty(rp)
  GR=[];
else
  GR=rp*gam.';
end
ch=(-gam.^2/2).^n;
cy=i*ch/pi;

nfac=factorial(n);
if nargout==2%% get coeff's of log(r) & poles:
  c_curlyY_log=2/nfac;
  sing_coeff(1)=c_curlyY_log*sum(A.*cy);

  c_curlyY_poles=-1/nfac;
  for j=0:n-1
    c_curlyY_poles=4*(j+(j==0))*(n-j)*gam.^(-2).*c_curlyY_poles;
    sing_coeff(j+2)=sum(A.*cy.*c_curlyY_poles);
  end
end

%% compute rest as a taylor series:
bn=-2/nfac/pi/i;
cn=-bn/2*( pi*i-2*log(1/2) + psi(1)+psi(n+1) );
curlyH0=cn+bn*log(gam);
Y=Y+sum(A.*ch.*curlyH0);%%limit as r->0 (OK m=0,1)

if isempty(GR)
  return
end

rr=1;
GR2=GR.^2;

Yb=0; Yc=0;
%GR2,max(max(GR2))
jj=find( GR2==max(max(GR2)) );
jj=jj(1);
j=0;
err=1;
while err>TOL
  j=j+1;
  rr=GR2.*rr;
  bn=-bn/4/j/(n+j);
  cn=-bn/2*( pi*i-2*log(1/2) + psi(j+1)+psi(n+j+1) );
  Yb=Yb+bn*rr;
  Yc=Yc+cn*rr;
  err=abs( cn*rr(jj)/Yc(jj) );
end

%  err_report={j,err};
%  if j==(MAXITS+1)
%    disp('warning [BES_hankel_series_smallr.m]: max no of iterations is reached'),MAXITS
%  end

Y(jp)=Y(jp) + ( Yc+(log(GR).*Yb) )*(A.*ch);

if 0
  plot( R,[real(Y0),imag(Y0)] );
  hold on, plot( [0 0],[real(Y0(1)),imag(Y0(1))],'.r' );
  xlim([-.1,max(r)+.1]),
  pause, close;
end
