function y=ND2(Z,hr,N,EE)
%% CALL: y=ND2(Z,hr,N);
%% Z={h,T,H_dim,th}, hr=h1/h0, N=no of imaginary roots req'd;
%% output={K,gam,H,th,[lam,-mu],char,L}.

h        = Z{1};
T        = Z{2};
H_dim    = Z{3};
thdeg    = Z{4};
H_infdep = 5;%% H>H_infdep is the infinite depth criterion

pram  = NDphyspram(0);
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
E        = EE(1,2);%% Pa

om    = 2*pi/T;
D     = E*h^3/12/(1-nu^2);
L_ice = (D/rho/g)^.25;
L5    = D/rho/om^2;
L     = L5^.2;

lam   = g/L/om^2;
mu    = rho_ice*h/rho/L;
del   = lam-mu;
H     = H_dim/L;
H     = min(H,H_infdep);

if iscell(N)
   decay_exp   = N{1};
   ndp         = N{2};
   tol         = 5*10^(-ndp-1);
   kap         = exp(-log(tol)/decay_exp);
   N           = round(kap/pi*H);
end

K     = RTS_ice_roots(del,H,N);
Er    = EE(1,1)/EE(1,2);
Dr    = Er*hr^3;
sigr  = hr*EE(2,1)/EE(2,2);

if sigr==1 & Dr==1
   gam = K;
else
   if sigr~=0
      L2    = Dr^.2;
      H2    = H/L2;
      del2  = (lam-sigr*mu)/L2;
      gam   = RTS_ice_roots(del2,H2,N)/L2;
   else
      gam   = RTS_wtr_roots(lam,H,N);
   end
end

char  = [L_ice sqrt(L_ice/g)];
y     = {K,gam,H,thdeg,[lam,-mu],char,L};
