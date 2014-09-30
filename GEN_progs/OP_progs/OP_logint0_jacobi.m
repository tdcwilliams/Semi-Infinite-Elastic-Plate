function Y=OP_logint0_jacobi(tt,alf,bet,M)
%% Y=\int_{-1}^t.w(s).log(t-s).ds,
%%  w(s)=(1-s)^\alpha.(1+s)^beta

TOL=1e-9;
Y=0*tt;
%% t==-1 => Y=0;

%% let u=log(t+1)-log(t-s)
%% => log(t-s)=log(t+1)-u;
%% => t-s=(t+1).e^{-u};
%% => s=t-(1+t).e^{-u};
%% => ds=(1+t).e^{-u}.du;

%% t==1 IS NEXT EASIEST:
j1=find(tt==1);
if ~isempty(j1)
  %% Y(1,alf,bet)
  %%  = \int_0^\infty.u^bet.e^{-(1+alf)*u}.ig1(u,alf,bet,M).du;
  %%  LET v=u/(1+alf) AND USE GAUSS-LAGUERRE QUADRATURE:
  Ngl=50;
  [vv,ww]=OP_numint_laguerre(bet,Ngl);
  gam=(1+M+alf);
  uu=vv/gam;
  ig=gam^(-1-bet)*ig1(uu,alf+M,bet);
  int=ww'*ig;
  err=10*TOL;
  while err>TOL
    int0=int;
    Ngl=Ngl+50;
    [vv,ww]=OP_numint_laguerre(bet,Ngl);
    uu=vv/gam;
    ig=gam^(-1-bet)*ig1(uu,alf+M,bet);
    int=ww'*ig;
    err=abs(1-int/int0);
  end
  Y(j1)=int;
end

%% NOW DO t\in(-1,1):
j0=find((tt>-1)&(tt<1));
if ~isempty(j0)
  Ngl=50;
  gam=(1+M);
  [vv,ww]=OP_numint_laguerre(bet,Ngl);
  uu=vv/gam;
  IG0=gam^(-1-bet)*ig0(tt(j0),uu,alf,bet);
  int1=IG0*ww;
  int2=IG0*(uu.*ww);
  err=10*TOL;
  while err>TOL
    int1_0=int1;
    int2_0=int2;
    Ngl=Ngl+50;
    [vv,ww]=OP_numint_laguerre(bet,Ngl);
    uu=vv/gam;
    IG0=gam^(-1-bet)*ig0(tt(j0),uu,alf,bet);
    int1=IG0*ww;
    int2=IG0*(uu.*ww);
    err=max( abs(1-int1/int1_0),abs(1-int2/int2_0) );
  end
  t=tt(j0);
  Y(j0)=(1+t).^(bet+1+M).*( log(1+t).*int1 - int2 );
end

function y=ig1(uu,alf,bet)
%% y(u,alf,bet)=2^(1+alf+bet)*(log(2)-u)*q(u)^bet
qq=1+0*uu;
jj=find(uu);
qq(jj)=( 1-exp(-uu(jj)) )./uu(jj);
y=2^(1+alf+bet)*(log(2)-uu).*qq.^bet;

function y=ig0(tt,uu,alf,bet)
%% y(t,u,alf,bet)=(1-t+e^{-u}*(1+t))^alf*q(u)^bet,
%%  q(u)=(1-e^{-u})/u;
[U,T]=meshgrid(uu,tt);
qq=1+0*U;
jj=find(U);
qq(jj)=( 1-exp(-U(jj)) )./U(jj);
y=( 1-T+exp(-U).*(1+T) ).^alf.*qq.^bet;