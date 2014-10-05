function T=GEN_get_ice_period(h,wavlen,H,E)
%% CALL: T=GEN_get_ice_period(h,wavlen,H,E)
%% calc's period corresponding to a given wavelength ('wavlen'),
%% given ice thickness 'h' and water depth (optional argument - if no
%% value is specified infinite water depth is assumed)

if h==0%% (ie if ice is water) use GEN_get_wtr_period.m:
   if ~exist('H')
      T  = GEN_get_wtr_period(wavlen);
   else
      T  = GEN_get_wtr_period(wavlen,H);
   end
   return;
end

pram  = NDphyspram(0);
if ~exist('E')
   E  = pram(1);%% Pa
end

g        = pram(2);%% m/s^2
rho      = pram(3);%% kg/m^3
rho_ice  = pram(4);%% kg/m^3
nu       = pram(5);

D     = E*h^3/12/(1-nu^2);
L_ice = (D/rho/g)^.25;
T_ice = sqrt(L_ice/g);

sig   = (rho_ice/rho)*(h/L_ice);
k     = 2*pi*L_ice/wavlen;
if nargin==2%%infinite depth
   kt = k;
else
   kt = k*tanh(k*H/L_ice);
end
x  = (k^4+1)*kt/(1+sig*kt);%% =(om*T_ice)^2
om = sqrt(x)/T_ice;
T  = 2*pi/om;
