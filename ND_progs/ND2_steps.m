function [Rts,HH,th,lam,sigsig,L]=ND2_steps(Z,hh,NN);
%% CALL: y=ND2(Z,hr,N);
%% Z={T,theta_inc,HH_dim}, hr=h1/h0, N=no of imaginary roots req'd;
%% output={K,gam,H,th,[lam,-mu],char,L}.

T=Z{1};
thdeg=Z{2};
th=thdeg*pi/180;
HH_dim=Z{3};
pram=NDphyspram(0);
H_infdep=5;%% H>H_infdep is the infinite depth criterion

E=pram(1);%% Pa
g=pram(2);%% m/s^2
rho=pram(3);%% kg/m^3
rho_ice=pram(4);%% kg/m^3
nu=pram(5);

hmax=max(hh);
hr=hh/hmax;
Dr=hr.^3;

om=2*pi/T;
D=E*hmax^3/12/(1-nu^2);
L_ice=(D/rho/g)^.25;
L=(D/rho/om^2)^.2;

lam=g/L/om^2;
sigsig=rho_ice*hh/rho/L;
HH=HH_dim/L;
%%
if 0%% apply infinite depth criterion:
  HH=H+max(sigsig)-sigsig;
  dH=max(HH)-H_infdep;
  HH=HH-dH*(dH>0);% [HH(2)-HH(1), sig]
end

%% GET ROOTS:
if hh(1)==hh(2) & HH(1)==HH(2)
  %% if both sides are the same,
  %% then only need to find one set of roots.
  N=max(NN);
  gam=RTS_ice_roots( lam-sigsig(1),HH(1),N );
  Rts={gam(1:NN(1)+3),gam(1:NN(2)+3)};
elseif prod(hh)==0
  %% then one side has open water:
  j0=find(hh==0);
  Rts{j0}=RTS_wtr_roots( lam,HH(j0),NN(j0) );

  %% and the other is ice-covered:
  j1=find(hh>0);
  Rts{j1}=RTS_ice_roots( lam-sigsig(j1),HH(j1),NN(j1) );
else
  %% both sides are different, but ice-covered:
  for j=1:2
    L_j=Dr(j)^.2;
    del_j=(lam-sigsig(j))/L_j;
    Rts{j}=RTS_ice_roots( del_j,HH(j)/L_j,NN(j) )/L_j;
  end
end