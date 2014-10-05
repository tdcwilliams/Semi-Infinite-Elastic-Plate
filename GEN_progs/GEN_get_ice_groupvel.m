function [ag,wavlen]=GEN_get_ice_groupvel(h,T,H,E)
%% CALL: [ag,wavlen]=GEN_get_ice_wavelength(h,T,varargin)
%% calc's group velocity (ag) & wavelength ('wavlen') corresponding to a given
%% ice thickness 'h', wave period 'T', and water depth
%% (optional argument - if no
%% value is specified infinite water depth is assumed)

DO_TEST  = 0;

if nargin==0
   DO_TEST  = 1;
   h        = 1;
   T        = 5;
   H        = 20;
   %%
   om    = 2*pi/T;
   eps   = 1e-6;
   omt   = om+eps;
   Tt    = 2*pi/omt;
elseif ~exist('H')
   H  = Inf;
end


pram  = NDphyspram(0);
if ~exist('E')
   E  = pram(1);%% Pa
end

g        = pram(2);%% m/s^2
rho      = pram(3);%% kg/m^3
rho_ice  = pram(4);%% kg/m^3
nu       = pram(5);
%%
om = 2*pi./T;

if h==0%% if ice is water, use GEN_get_wtr_period.m
   if H==Inf
      k  = om.^2/g;
      ap = om./k;
      %%
      ag       = ap/2;
      wavlen   = 2*pi./k;
   else
      wavlen   = GEN_get_wtr_wavelength(T,H);
      k        = 2*pi./lam;
      %%
      tk = tanh(k*H);
      sk = sech(k*H).^2;
      ag = g./(2*om).*(tk+H*k.*sk);
   end
   return;
end

D     = E*h.^3/12/(1-nu^2);
L     = ((D/rho)./om.^2).^.2;
%  char=[L,L_ice];
lam   = g./L./om.^2;
mu    = (rho_ice*h/rho)./L;
del   = lam-mu;

if H==Inf
   guess = L.*om.^2/g;
   k     = gen_root_ice(del,Inf,guess)./L; 
   ag    = ( 5*D*k.^4 + rho*g-rho_ice*h*om.^2 )./...
            ( 2*om.*(rho+rho_ice*k*h) );
else
   guess = L.*om.^2/g;
   k     = gen_root_ice(del,H./L,guess)./L; 
   %%
   kH    = k.*H;
   tk    = tanh(kH);
   sk    = sech(kH).^2;
   Lam   = D*k.^4+rho*g-rho_ice*h*om.^2;
   %%
   ag = ( (4*D*k.^4+Lam).*tk+Lam.*kH.*sk )./...
            ( 2*om.*(rho+rho_ice*h*k.*tk) );
end

wavlen   = 2*pi./k;

if DO_TEST
   lam   = GEN_get_ice_wavelength(h,T,H,E);
   [lam,wavlen]
   lamt  = GEN_get_ice_wavelength(h,Tt,H,E);
   k  = 2*pi/lam;
   kt = 2*pi/lamt;
   [ag,eps/(kt-k)]
end


function k=gen_root_ice(del,H,guess)
%% finds the root of the ice dispersion relation nearest to 'guess'.

tol   = 1e-8;
k0    = guess;
dk    = NR_corr_term(k0,del,H);
k     = k0-dk;
while max(abs(dk)) > tol
   k0 = k;
   dk = NR_corr_term(k0,del,H);
   k  = k0-dk;
end

function dk=NR_corr_term(k,del,H)
%% dk=f/f_k, where f has the same zeros as of the dispersion function, 
%% is the correction term in the Newton-Rhapson method for finding zeros in f.

Lam   = k.^4+del;
Lampr = 5*k.^4+del;
x     = 7.5;
if H==Inf
   f  = Lam.*k-1;
   df = Lampr;
else
   kH = k.*H;
   f  = Lam.*k.*sinh(kH)-cosh(kH);
   df = Lam.*kH.*cosh(kH)+(Lampr-H).*sinh(kH);
   %%
   jj = find(real(kH)>x);

   if ~isempty(jj)
      f(jj)    = Lam(jj).*k(jj).*tanh(kH(jj))-1;
      df(jj)   = Lam(jj).*kH(jj)+...
                  (Lampr(jj)-H(jj)).*tanh(kH(jj));
   end
end
dk = f./df;
