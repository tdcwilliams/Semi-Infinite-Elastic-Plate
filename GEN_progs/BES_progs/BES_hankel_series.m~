function [Y,sing_coeff]=BES_hankel_series(r,A,gam,n,keepsing)
%% call: [Y,log_coeff]=BES_hankel_series(x,An,gam,m,keepsing)
%% calc's Y=Lw=(1/r.d/dr)^nW, where
%% W=sum_{n=-2}^length(gam).A_n.H_0(gam_n|x|)
%% uses the result that L.H_0(kr)=(-k/r)^nH_n(kr)
%% used to calc the thin plate 'in vacuo' Green's function
%% or the full floating thin plate Green's function at the surface.

sing_coeff=zeros(n+1,1);
sz=size(r);
nr=prod(sz);
R=zeros(nr,1);%% want to work with a column vector;
Y0=R;
R(1:nr)=r(1:nr);
Y=zeros(sz);

if nargin==4%% default is subtract off singularities
  keepsing=0;
elseif keepsing==1%% don't subtract them off.
  GR=R*gam.';
  Y0=(1./GR).^n.*besselh(n,1,GR);
  C=A.*(i*gam).^(2*n);
  Y(1:nr)=Y0*C;
  return
end

%% Calc values for small r using a power series:
jm=find(R<=1);
Rm=R(jm);
if isempty(Rm)==0
  %{Rm,isempty(Rm)}
  Y0(jm)=BES_hankel_series_smallr(Rm,A,gam,n,keepsing);
end

%% Calc values for larger r by calculating hankel functions
%% directly and subtracting the singularities:
jp=find(R>=0.5);
Rp=R(jp);
if ~isempty(Rp)
  Ysing=zeros(size(Rp));

  ch=(-gam.^2/2).^n;
  cy=i*ch/pi;

  nfac=factorial(n);
  c_curlyY_log=2/nfac;
  sing_coeff(1)=c_curlyY_log*sum(A.*cy);
  Ysing=sing_coeff(1)*log(Rp);

  c_curlyY_poles=-1/nfac;
  rr=1;
  for j=0:n-1
    rr=rr./Rp.^2;
    c_curlyY_poles=...
      4*(j+(j==0))*(n-j)*gam.^(-2).*c_curlyY_poles;
    sing_coeff(j+2)=sum(A.*cy.*c_curlyY_poles);
    Ysing=Ysing+sing_coeff(j+2)*rr;
  end
  GR=Rp*gam.';
  Y0(jp)=GR.^(-n).*besselh(n,1,GR)*(A.*(i*gam).^(2*n)) - Ysing;
end
Y(1:nr)=Y0;

if 0
  plot( R,[real(Y0),imag(Y0)] );
  hold on, plot( [0 0],[real(Y0(1)),imag(Y0(1))],'.r' );
  xlim([-.1,max(r)+.1]),
  pause, close;
end