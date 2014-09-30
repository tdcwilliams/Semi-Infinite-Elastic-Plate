function [L,lam,gam0,alpy]=NDbasics_fluid(Z)
%% CALL: L=NDbasics(Z,hr)
%% del0=[lam,-mu]
%% INPUT:
%% Z=T or [T th] or [T th L]

T=Z(1);
th=0;
n=length(Z);
wantL=(n<3);
if n==2
  th=pi*Z(2)/180;
elseif n==3
  th=pi*Z(2)/180;
  L=Z(3);
end
%%

g=9.81;
om=2*pi/T;
gam0_dim=om^2/g;
if want_L
  L=1/gam0_dim;
end

gam0=gam0_dim*L;
lam=1/gam0;
alpy=gam0*sin(th);