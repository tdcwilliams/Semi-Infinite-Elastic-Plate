function [del0, L, xtra]=NDmostbasic(Z,hr)
%% CALL: [del0, L, gam, char]=NDmostbasic(Z,hr)
%% del0=[lam,-mu], xtra=[th,nu]
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
L=(D/rho/om^2)^.2;
lam=g/L/om^2; mu=rho_ice*h/rho/L;
del0=[lam,-mu];
xtra=[th,nu];
