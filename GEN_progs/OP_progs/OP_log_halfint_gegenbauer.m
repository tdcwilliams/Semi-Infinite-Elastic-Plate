function Y=OP_logint1_gegenbauer_even(tt,alf,Nterms)
%% CALL: Y=OP_logint1_gegenbauer_even(tt,alf,Nterms)
%% Y=2*\int_0^1.(1-s^2)^(\alf-.5).C_{2n}^(\alf)(s).log|t-s|ds

if 0
    F2=inline('( 2*alpha*(1+alpha)*((.5+sigma/2).^2-(.5+tau/2).^2)-4*alpha*(1+alpha)*(.5+tau/2).*(sigma-tau)/2 ).*(1-sigma).^(alpha-.5).*(.75+sigma/4).^(alpha-.5).*log(abs(tau/2-sigma/2))');

F2a=inline('(1-sigma).^(alpha-.5).*(.75+sigma/4).^(alpha-.5).*log(abs(tau/2-sigma/2))')
F2b=inline('2*alpha*(1+alpha)*((.5+sigma/2).^2-(.5+tau/2).^2)-4*alpha*(1+alpha)*(.5+tau/2).*(sigma-tau)/2')

  tau=.1; t=.5+tau/2;
  sig=.2; s=.5+sig/2;
  alf2=alf-.5;
  pj_t=2*alf*(alf+1)*t^2-alf;
  dpj_t=4*alf*(alf+1)*t;
  pj_s=2*alf*(alf+1)*s^2-alf;
  tst_ig=[F2(alf,sig,tau), (1-sig).^alf2.*ig_n(sig,tau,alf2,pj_t,dpj_t,pj_s)]
  CC=pj_s-pj_t-dpj_t*(s-t);
  [F2b(alf,sig,tau),CC]
  [CC*(1-s^2)^alf2*log(abs(t-s)),CC*(1-s^2)^alf2*log(abs(tau-sig)/2), CC*(1-sig)^alf2*(.75+sig/4)^alf2*log(abs(tau-sig)/2)]
  return
end

TOL=1e-9;
Y=zeros(length(tt),Nterms+1);

%% 1ST CALC n=0:
tau=2*tt-1;
alf2=alf-.5;
%% I0=\int_{-1}^1\log|\tau-\sig|*(1-\sig)^(\alf-.5)d\sig
M=0;
I0=OP_logint0_jacobi(tau,alf2,0,M)+...
      OP_logint0_jacobi(-tau,0,alf2,M);
Y(:,1)=( (3+tau)/4 ).^alf2.*I0;

%% I1=\int_{-1}^1\log|\tau-\sig|*(\tau-\sig)*...
%%              (1-\sig)^(\alf-.5)d\sig
M=1;
I1=OP_logint0_jacobi(tau,alf2,0,M)-...
      OP_logint0_jacobi(-tau,0,alf2,M);
Y(:,1)=Y(:,1)-.25*alf2*( (3+tau)/4 ).^(alf2-1).*I1;

%% 'HARD' INTEGRAL FOR n=0 :)
Nint=0;
int=1e6;
err=10*TOL;
while err>TOL
  int0=int;
  Nint=Nint+50;
  [sig,wJ]=OP_numint_jacobi(alf2,0,Nint);
  int1=-log(2)*sum( wJ.*( (3+sig)/4 ).^alf2 );
  %%
  int = int1 + ig0(sig,tau,alf2)*wJ;
  err=max(abs(1-int0/int));
end
Y(:,1)=Y(:,1)+int;

if Nterms>=1%% do n=1;
  Y(:,2)=2*alf*tt.*Y(:,1) - alf*tau*int1 +...
    - alf*((3+tau)/4).^alf2.*I1;
  Nint=0;
  int=1e6;
  err=10*TOL;
  while err>TOL
    int0=int;
    Nint=Nint+50;
    [sig,wJ]=OP_numint_jacobi(alf2,0,Nint);
    int2=-alf*log(2)*sum( wJ.*sig.*( (3+sig)/4 ).^alf2 );
    %%
    int = int2 + alf*ig1(sig,tau,alf2)*wJ;
    err=max(abs(1-int0/int));
  end
  Y(:,2)=Y(:,2)+int;
end

%% INTEGRALS INVOLVING ODD POLYNOMIALS CONVERGE FASTER
%%  THAN FOR EVEN ONES => IF THE HIGHEST POLYNOMIAL IS ODD
%%   & YOU JUST USE ENOUGH POINTS TO MAKE IT CONVERGE,
%%    IT MAY NOT BE ENOUGH FOR THE LOWER EVEN INTEGRALS:
WNT_XTRA=( ~( round(Nterms/2)==Nterms/2 ) & Nterms>=2);
if WNT_XTRA
  Nterms=Nterms+1;
end

if Nterms>=2
  for j=Nterms%% GET RESULTS TO CONVERGE FOR HIGHEST (EVEN)
              %%  POLYNOMIAL FIRST
              %%   => DON'T NEED 'while' FOR LOWER ONES:
    Cn_vals_t=OP_interp_gegenbauer(tt,alf,{Nterms});
    dCn_vals_t=2*alf*OP_interp_gegenbauer(tt,alf+1,{Nterms-1});
    %%
    pj_t=Cn_vals_t(:,j+1);
    dpj_t=dCn_vals_t(:,j);
    Y(:,j+1)=pj_t.*Y(:,1) +dpj_t.*( Y(:,2)/2/alf-tt.*Y(:,1) );
    %%
    Nint=0;
    int=1e6;
    err=10*TOL;
    while err>TOL
      int0=int;
      Nint=Nint+50;
      [sig,wJ]=OP_numint_jacobi(alf2,0,Nint);
      Cn_vals_s=OP_interp_gegenbauer((1+sig)/2,alf,{Nterms});
      pj_s=Cn_vals_s(:,j+1);
      %%
      int = ig_n(sig,tau,alf2,pj_t,dpj_t,pj_s)*wJ;
      err=max(abs(1-int0/int));
    end
    Y(:,j+1)=Y(:,j+1)+int;
  end
  %%
  for j=2:Nterms-1
    pj_t=Cn_vals_t(:,j+1);
    dpj_t=dCn_vals_t(:,j);
    Y(:,j+1)=pj_t.*Y(:,1) +dpj_t.*( Y(:,2)/2/alf-tt.*Y(:,1) );
    %%
    pj_s=Cn_vals_s(:,j+1);
    int = ig_n(sig,tau,alf2,pj_t,dpj_t,pj_s)*wJ;
    Y(:,j+1)=Y(:,j+1)+int;
  end
end

if WNT_XTRA
  Y(:,Nterms+1)=[];
end

if 1
  F=inline('2*(1-s.^2).^(alpha-.5).*log(abs(t-s))');
  II0=quad(@(s)F(alf,s,tt(1)),0,1,1e-11)
  %%
  F=inline('2*(2*alpha*s).*(1-s.^2).^(alpha-.5).*log(abs(t-s))');
  II1=quad(@(s)F(alf,s,tt(1)),0,1,1e-11)
  %%
  F=inline('2*(2*alpha*(1+alpha)*s.^2-alpha).*(1-s.^2).^(alpha-.5).*log(abs(t-s))');
  II2=quad(@(s)F(alf,s,tt(1)),0,1,1e-11)
  %%
  F2=inline('( 2*alpha*(1+alpha)*((.5+sigma/2).^2-(.5+tau/2).^2)-4*alpha*(1+alpha)*(.5+tau/2).*(sigma-tau)/2 ).*(1-sigma).^(alpha-.5).*(.75+sigma/4).^(alpha-.5).*log(abs(tau/2-sigma/2))');
  %%
  F3=inline('2*( 2*alpha*(1+alpha)*(s.^2-t.^2)-4*alpha*(1+alpha)*t.*(s-t) ).*(1-s.^2).^(alpha-.5).*log(abs(t-s))');

%   INT=[quad(@(s)F3(alf,s,tt(1)),0,1,1e-11),quad(@(sig)F2(alf,sig,tau(1)),-1,1,1e-11),int]

%    INT+pj_t.*Y(:,1) +dpj_t.*( Y(:,2)/2/alf-tt.*Y(:,1) )
  %%
%    tst_ig=[F2(alf,sig(1),tau(1)), (1-sig(1)).^alf2.*ig_n(sig(1),tau(1),alf2,pj_t(1),dpj_t(1),pj_s(1))]

%    [pj_t(1),2*alf*(alf+1)*tt(1)^2-alf]
%    [pj_s(1),2*alf*(alf+1)*(.5+sig(1)/2)^2-alf]
%    [dpj_t(1),4*alf*(alf+1)*tt(1)]


end
return;
%  F1=inline('(1-sig).^(alpha-.5).*( (3+sig)/4 ).^(alpha-.5).*log(.5)')
%  tst_int1=[int1,quad(@(sig)F1(alf,sig),-1,1,1e-8)]

%  F2=inline('(1-s).^(alpha-.5).*log(abs(t-s)).*( ((3+s)/4).^(alpha-.5) - ((3+t)/4).^(alpha-.5)-.25*(alpha-.5)*(s-t).*((3+t)/4).^(alpha-1.5) )')
%  IF2=[quad(@(sig)F2(alf,sig,tau(1)),-1,1,1e-8), int2]
%  if2b=[quad(@(sig)ig(sig,tau(1),alf2),-1,1,1e-8),int2]
%%


%IF2B=[quad(@(sig)ig(sig,tau(1),alf2),-1,1,1e-8)+int1, int]



%%
%  
%  G=inline('(1-sig).^(alpha-.5).*( (3+sig)/4 ).^(alpha-.5).*log(abs(tau-sig)/2)')
%  INT=quad(@(sig)G(alf,sig,tau(1)),-1,1,1e-10)
%  INTb=[IF2(1)+int1(1)+((3+tau(1))/4)^alf2*I0(1)-.25*alf2*((3+tau(1))/4)^(alf2-1)*I1(1),...
%  IF2(1)+int1(1)+INT0+INT1B]

%  [,((3+tau(1))/4)^alf2*I0(1)-4*alf2*((3+tau(1))/4)^(alf2-1)*I1(1)]












tau=2*tt-1;
[C2n_vals_t,h2nC]=OP_interp_gegenbauer(tt,alf,{2*Nterms});
jj=2:2:(2*Nterms+1);
h2nC(jj)=[];
C2n_vals_t(:,jj)=[];
alf2=alf-.5;
eN=find( (0:2*Nterms)'==2*Nterms );
Y=0*C2n_vals_t;
%%
for j=Nterms
  Nint=50;
  [sig,wJ]=OP_numint_jacobi(alf-.5,0,Nint);
  ss=.5*(1+sig);
  C2n_vals_s=OP_interp_gegenbauer(ss,alf,eN);
  %%
  [Sig,Tau]=meshgrid(sig,tau);
  [S,T]=meshgrid(ss,tt);
  IG=0*Tau;
  jp=find(Sig~=Tau);
  %%
  Ds=diag( C2n_vals_s/h2nC(j+1) );
  Dt=diag( C2n_vals_t(:,j+1)/h2nC(j+1) );
  IG0=( (3+Sig)/4 ).^alf2*Ds-Dt*( (3+Tau)/4 ).^alf2;
  IG(jp)=IG0(jp).*log( abs(T(jp)-S(jp))/2 );
  int1=IG*wJ;
  ERR=10*TOL;
  while ERR>TOL
    Nint=Nint+50
    int10=int1;
    %%
    [sig,wJ]=OP_numint_jacobi(alf-.5,0,Nint);
    ss=.5*(1+sig);
    C2n_vals_s=OP_interp_gegenbauer(ss,alf,eN);
    %%
    [Sig,Tau]=meshgrid(sig,tau);
    [S,T]=meshgrid(ss,tt);
    IG=0*Tau;
    jp=find(Sig~=Tau);
    %%
    Ds=diag( C2n_vals_s/h2nC(j+1) );
    IG0=( (3+Sig)/4 ).^alf2*Ds-Dt*( (3+Tau)/4 ).^alf2;
    IG(jp)=IG0(jp).*log( abs(T(jp)-S(jp))/2 );
    int1=IG*wJ
    %%
    err=abs(1-int10./int1)
    ERR=max(err);
  end
end

int2=OP_logint0_jacobi(tau,alf-.5,0)+...
      OP_logint0_jacobi(-tau,0,alf-.5)
        -log(2)*sum(wJ);
Y(:,j+1)=int1+...
    ((3+tau)/4).^alf2.*C2n_vals_t(:,j+1)/h2nC(j+1).*int2;


%% ASSUME THAT INTEGRALS FOR SMALLER j WILL CONVERGE
%%  MORE QUICKLY THAN j=Nterms (SO NO NEED TO REPEAT
%%  'while' LOOP):
C2n_vals_s=OP_interp_gegenbauer(ss,alf,{2*Nterms+1});
C2n_vals_s(:,jj)=[];
for j=0:Nterms-1
  Ds=diag( C2n_vals_s(:,j+1)/h2nC(j+1) );
  Dt=diag( C2n_vals_t(:,j+1)/h2nC(j+1) );
  IG0=( (3+Sig)/4 ).^alf2*Ds-Dt*( (3+Tau)/4 ).^alf2;
  IG(jp)=IG0(jp).*log( abs(T(jp)-S(jp))/2 );
  int1=IG*wJ;
  %%
  Y(:,j+1)=int1+...
    ((3+tau)/4).^alf2.*C2n_vals_t(:,j+1)/h2nC(j+1).*int2;
end

function Y=ig0(sig,tau,alf2)

[Sig,Tau]=meshgrid(sig,tau);
Y=0*Sig;
jnz=find(Sig~=Tau);
%%
y0=(3+Sig).^alf2 - (3+Tau).^alf2 +...
     -alf2*(Sig-Tau).*(3+Tau).^(alf2-1);
Y(jnz)=y0(jnz).*log(abs(Tau(jnz)-Sig(jnz)))/4^alf2;

function Y=ig1(sig,tau,alf2)

[Sig,Tau]=meshgrid(sig,tau);
Y=0*Sig;
jnz=find(Sig~=Tau);
%%
y0=( (3+Sig).^alf2 - (3+Tau).^alf2 ).*(Sig-Tau);
Y(jnz)=y0(jnz).*log(abs(Tau(jnz)-Sig(jnz)))/4^alf2;

function Y=ig_n(sig,tau,alf2,pj_t,dpj_t,pj_s)

[Sig,Tau]=meshgrid(sig,tau);
[pj_S,pj_T]=meshgrid(pj_s,pj_t);
Y=0*Sig;
jnz=find(Sig~=Tau);
%%
y0=pj_S-pj_T-diag(dpj_t/2)*(Sig-Tau);
Y(jnz)=(.75+Sig(jnz)/4).^alf2.*...
         log(abs(Tau(jnz)-Sig(jnz))/2).*...
           y0(jnz);