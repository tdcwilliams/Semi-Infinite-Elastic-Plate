function period=NDget_periodfromlam(h,lam)

%%
pram=NDphyspram(0);
E=pram(1); %Pa
g=pram(2); %m/s^2
rho=pram(3); %kg/m^3
rho_ice=pram(4); %kg/m^3
nu=pram(5);
%%
bet=12*(1-nu^2)*rho/E/h^3;
om=(bet*g^5)^(1/8)*lam.^(-5/8);
period=2*pi./om;