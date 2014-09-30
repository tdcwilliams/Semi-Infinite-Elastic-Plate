function [del0, L, gam, el, char,pram]=NDbasics(Z,hr)
%% CALL: [del0, L, gam, char]=NDbasics(Z,hr)
%% del0=[lam,-mu]
%% INPUT:
%% Z=T or [T th] or [T th h]

T=Z(1);
th=0;
h=1;
n=length(Z);
if n==2
  th=pi*Z(2)/180;
elseif n==3
  th=pi*Z(2)/180;
  h=Z(3);
end
%%
if nargin==1
  hr=1;
end
%%
pram=NDphyspram(0);
E=pram(1); %Pa
g=pram(2); %m/s^2
rho=pram(3); %kg/m^3
rho_ice=pram(4); %kg/m^3
nu=pram(5);
%%
om=2*pi/T;
D=E*h^3/12/(1-nu^2);
L_ice=(D/rho/g)^.25;
%
L=(D/rho/om^2)^.2%{L,(E/rho/12/(1-nu^2))^.2*h^.6*om^(-.4)}
char=[L_ice sqrt(L_ice/g)];
lam=g/L/om^2; mu=rho_ice*h/rho/L;
del0=[lam,-mu];

if 1
  del=sum(del0)
  alp=(E/rho/12/(1-nu^2))^.2;
  {alp*del*h^.6*om^1.6,g-.9*h*om^2}
end


if nargout>=3
  gam=zeros(5,1);
  r=roots([hr^3 0 0 0 lam-hr*mu -1]);
  gam(1)=r( find(r>0 & imag(r)==0) ); el=gam(1)*sin(th);
  gam(2:3)=r( find(r>0 & imag(r)~=0) );
  gam(4:5)=r (find(r<0) );
end
