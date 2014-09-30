function [period,L]=NDreverse_ice_nondim_tp(del,h)

%%
pram=NDphyspram(0);
E=pram(1); %Pa
g=pram(2); %m/s^2
rho_wtr=pram(3); %kg/m^3
rho_ice=pram(4); %kg/m^3
nu=pram(5);
%%
alp=(E/rho_wtr/12/(1-nu^2))^.2;
rho=rho_ice/rho_wtr;
rts=roots([rho*h,alp*del*h^.6,0,0,0,-g]);
om_to_pt4=rts(find(imag(rts)==0 & rts>0));
om=om_to_pt4^2.5;
%%
period=2*pi/om;
L=alp*h^.6/om^.4;