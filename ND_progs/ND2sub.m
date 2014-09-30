function [Rts,HHsig,th,del0,L] = ND2sub(Z,hh,NNN,EE,rho);
%% CALL: [Rts,HHsig,th,del0,L] = ND2sub(Z,hh,NNN,EE,rho);
%% Z  = {T,H_dim,th}
%%    T     = period;
%%    H_dim = depth (relative to the most submerged ice);
%%    th    = angle of incidence;
%% hh = [h1,h2] 
%%  - ice thickness on left/right;
%% NNN   = no of imaginary roots req'd
%%  - N1 or [N1,N2], or [N1,N2,N3];
%% EE = E1 or [E1,E2] or [E1,E2;rho1,rho2;nu1,nu2];
%%    Young's modulus, ice density, Poisson's ratio;
%% rho is the water density;

DO_TEST  = 0;
if nargin==0%%test inputs
   period   = 10;
   thdeg    = 0;
   H_dim    = 100;
   Z        = {period,thdeg,H_dim};
   %%
   hh = [1 2];
   %hh = [0 1];
   %%
   NNN   = [10 10];
   %%
   E0 = 1*NDphyspram(1);
   if 1
      EE = [1*E0,2*E0;
            922.5,922.5;
            .3,.3];
   else
      EE = E0;
   end
   %%
   rho      = 1025;
   DO_TEST  = 1;
end

T        = Z{1};
thdeg    = Z{2};
th       = thdeg*pi/180;
H_dim    = Z{3};
H_infdep = 5;%% H>H_infdep is the infinite depth criterion

pram  = NDphyspram(0);%[E,g,rho,rho_ice,nu]
g     = pram(2);%% m/s^2

if ~exist('rho')
   rho   = pram(3);%% kg/m^3
end

if ~exist('EE')
   E        = pram(1);%% Pa
   rho_ice  = pram(4);%% kg/m^3
   nu       = pram(5);
   EE       = [E,E;rho_ice,rho_ice;nu,nu];
elseif prod(size(EE))==1
   %%if only a scalar is entered, this is Young's modulus for any ice
   E        = EE(1);
   rho_ice  = pram(4);
   nu       = pram(5);
   EE       = [E,E;rho_ice,rho_ice;nu,nu];
end

rho_ice  = pram(4);%% kg/m^3
nu       = pram(5);

rhorho   = EE(2,:);
nunu     = EE(3,:);
EE       = EE(1,:);

h1 = hh(1);
E1 = EE(1);
h2 = hh(2);
E2 = EE(2);
D1 = E1*h1^3/12/(1-nu^2);
D2 = E2*h2^3/12/(1-nu^2);
%%
rho1  = rhorho(1);
rho2  = rhorho(2);
sig1  = rho1*h1/rho;
sig2  = rho2*h2/rho;

if sig1>sig2%%when submergence is involved,
        %%geometry (h) is most impt, E less so
   D     = D1;%%ref D (max D)
   jmax  = 1;
   jmin  = 2;
else%%default: rhs is ref
   D     = D2;
   jmax  = 2;
   jmin  = 1;
end
Dr = [D1 D2]/D;

om    =  2*pi/T;
L_ice = (D/rho/g)^.25;
L     = (D/rho/om^2)^.2;

lam      = g/L/om^2;
sigsig   = (rhorho.*hh)/rho/L;
H        = H_dim/L;
%%
sig   = abs(sigsig*[1;-1]);%%difference in bottom levels;
HH    = H+max(sigsig)-sigsig;
dH    = max(HH)-H_infdep;
HH    = HH-dH*(dH>0);% [HH(2)-HH(1), sig]
HHsig = [HH,sig];
mu    = max(sigsig);
del0  = [lam,-mu];
%%
if iscell(NNN)
   decay_exp   = NNN{1};
   ndp         = NNN{2};
   tol         = 5*10^(-ndp-1);
   kap         = exp(-log(tol)/decay_exp);
   Nreq        = round(kap/pi*HHsig);
   N0          = Nreq(1);
   N1          = Nreq(2);
   N_sig       = Nreq(3);
   WANT3       = (length(NNN)==4);
   NNN         = Nreq;
else
   N0    = NNN(1);
   N1    = NNN(2);
   WANT3 = (length(NNN)==3);
   if WANT3
      N_sig = NNN(3);
   end
end

if hh(1)==hh(2) & D1==D2
   gam   = RTS_ice_roots(lam-mu,HH(1),N0);
   Rts   = {gam,gam,[]};
else
   %jmax  = find(hh==hmax);
   %jmin  = find(hh~=hmax);
   %[lam-mu,HH(jmax),NNN(jmax)]
   Rts{jmax}   = RTS_ice_roots( lam-mu,HH(jmax),NNN(jmax) );%Rts{jmax},pause
   if Dr(jmin)==0
      Rts{jmin}   = RTS_wtr_roots( lam,HH(jmin),NNN(jmin) );
      if WANT3
         Rts{3}   = RTS_wtr_roots(lam,sig,N_sig);
      end
   else
      L_j         = Dr(jmin)^.2;
      del_j       = (lam-sigsig(jmin))/L_j;
      Rts{jmin}   = RTS_ice_roots( del_j,HH(jmin)/L_j,NNN(jmin) )/L_j;
      if WANT3
         Rts{3}   = RTS_ice_roots(del_j,sig/L_j,N_sig)/L_j;
      end
   end
end

if DO_TEST
   k1 = Rts{1}/L
   k2 = Rts{2}/L
end
