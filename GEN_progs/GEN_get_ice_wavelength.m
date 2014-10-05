function [wavlen,char]  = GEN_get_ice_wavelength(h,T,H,E)
%% CALL: wavlen=GEN_get_ice_wavelength(h,T,varargin)
%% calc's wavelength ('wavlen') corresponding to a given
%% ice thickness 'h', wave period 'T', and water depth
%% (optional argument - if no
%% value is specified infinite water depth is assumed)

do_test  = 0;
if nargin==0
   do_test  = 1;
   h        = 1;
   T        = 5;
   H        = 100;
   E        = NDphyspram(1);
end

if h==0%% if ice is water, use GEN_get_wtr_period.m
   if ~exist('H');
      wavlen   = GEN_get_wtr_wavelength(T);
   else
      wavlen   = GEN_get_wtr_wavelength(T,H);
   end
   char  = [];
   return;
end

if ~exist('H')
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

om       = 2*pi./T;
D        = E*h.^3/12/(1-nu^2);
L        = ((D/rho)./om.^2).^.2;
lam      = g./L./om.^2;
mu       = (rho_ice*h/rho)./L;
del      = lam-mu;
wavlen   = 0*del;
%%
L_ice =(D/rho/g).^.25;
char  = {L,L_ice};


%% infinite depth results:
for j=1:length(del)
   r  = roots([1 0 0 0 del(j) -1]);
   k  = r(find(r>0 & imag(r)==0));
   if H~=Inf
      %% if want wavelength for finite depth,
      %% use the infinite depth 'k' as a seed.
      k  = gen_root_ice(del(j),H/L(j),k);
   end
   wavlen(j)   = 2*pi*L(j)/k;
end

if do_test
   disp('test dispersion relation for T, h, E:');
   disp({T,h,E});
   %%
   kk       = linspace(0,2*k,51);
   [dk,ff]  = NR_corr_term(kk,del,H/L);
   [dk,f0]  = NR_corr_term(k,del,H/L);
   [f0,disprel(k,1,del,H/L)*cosh(k*H/L),disprel(k/L,L^5,del*L,H)*cosh(k*H/L)]
   %f=(k^4+0.424104485497080)*k*tanh(k*7.98351132886686)-1
   %%
   %k,1,del,H/L,L
   %f0 = disprel(k,1,del,H/L)
   %%
   plot(kk,ff,'b',...
        kk,0*kk,'k',...
        k,0,'.r');
   ylim([-1 1]*5);
   xlim([0 2]*k);
end

function k=gen_root_ice(del,H,guess)
%% finds the root of the ice dispersion relation nearest to 'guess'.

tol   = 1e-12;
k0    = guess;
dk    = NR_corr_term(k0,del,H);
k     = k0-dk;
while abs(dk) > tol
   k0 = k;
   dk = NR_corr_term(k0,del,H);
   k  = k0-dk;
end

function [dk,f] = NR_corr_term(k,del,H)
%% dk=f/f_k, where f has the same zeros as of the dispersion function, 
%% is the correction term in the Newton-Rhapson method for finding zeros in f.

Lam   = k.^4+del;
Lampr = 5*k.^4+del;
x     = 7.5;
if real(k*H)<=x
   f  = Lam.*k.*sinh(k*H)-cosh(k*H);
   df = Lam.*(k*H).*cosh(k*H)+(Lampr-H).*sinh(k*H);
else
   f  = Lam.*k.*tanh(k*H)-1;
   df = Lam.*k*H+(Lampr-H).*tanh(k*H);
end
dk = f./df;

function [f,df]   = disprel(gam,c5,c1,H)
kt = gam.*tanh(gam*H);
%%
f  = (c5*gam.^4+c1).*kt-1;
df = (5*c5*gam.^4+c1).*tanh(gam*H)+(c5*gam.^4+c1).*(gam*H).*sech(gam*H).^2;
