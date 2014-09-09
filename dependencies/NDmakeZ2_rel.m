function [Rts,rts,H,el,del0,L,char]=...
            NDmakeZ2_rel(Z,hh,N,EE,rho)
%% MAIN NONDIMENSIONALISATION PROGRAM WHEN 2 SHEET THICKNESSES ARE INVOLVED.
%% Differs from NDmakeZ2.m in that it nondimensionalises with respect to
%% the thickest ice sheet,
%% whereas the other always ND's wrt to the lh ice sheet.
%%
%% CALL: [Rts,rts,H,el,del0,L,char]=NDmakeZ2_rel(Z,[h0 h1],N)
%% INPUTS: Z=T or [T,th] or {T,th,H_dim},
%% hr=h1/h0 (hr=0 => open water to far right),
%% N is the number of imaginary roots wanted.
%%
%% OUTPUTS: Rts={K,gam}, where K is the (N+3)-vector of the roots
%% satisfying the dispersion relation for the lhs, gam is the
%% (N+3)- or )N+1)-vector satisfying the dispersion relation for the rhs;
%% el=K(1)*sin(theta);
%% rts={kn,alp}, where kn=sqrt(K.^2-el^2) and alp=sqrt(gam.^2-el^2)
%% are the wavenumbers in the x direction
%% (square roots are taken from the upper complex half-plane);
%% H is the nondimensional depth, del0=[lam,-mu], L is the natural length
%% and char=[L_ice, T_ice], where L_ice is the characteristic length
%% for the lh ice and T_ice is its characteristic time.

do_test  = 0;
if nargin==0%%do a test
   do_test  = 1;
   hh       = [1 2];
   %youngs   = [5e9 6e9];
   %youngs   = [5e9 5e9];
   youngs   = NDphyspram(1)*[1 .8];
   N        = 10;
   %%
   T     = 5;
   theta = 0;
   H_dim = 100;
   Z     = {T,theta,H_dim};
end

if length(Z)==1
   T     = Z;
   Z     = NDphyspram([]);
   Z{2}  = T;
elseif length(Z)==2
   T     = Z(1);
   theta = Z(2);
   Z     = NDphyspram([]);
   Z{2}  = T;
   Z{4}  = theta;
elseif length(Z)==3
   T     = Z{1};
   theta = Z{2};
   H_dim = Z{3};
   Z     = {1,T,H_dim,theta};
end

pram  = NDphyspram(0);%[E,g,rho,rho_ice,nu]
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
if ~exist('rho')
   rho   = 1025;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get impt parameters:
[hh2,jj] = sort(hh);
hmax     = hh(jj(2));
Z{1}     = hmax;%testZ=Z
hr       = hh(jj(1))/hmax;
EE       = youngs(:,jj);
y        = ND2(Z,hr,N,EE,rho);%% y={K,gam,H,th,del0,char,L}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H     = y{3};
th    = pi*y{4}/180;
del0  = y{5};
char  = y{6};
L     = y{7};

%% sort out roots:
Rts{jj(2)}  = y{1};%%larger thickness
Rts{jj(1)}  = y{2};%%smaller thickness
gam1        = Rts{1};
gam2        = Rts{2};
el          = gam1(1)*sin(th);
%%
alp1  = sqrt(gam1.^2-el^2);
alp1  = sign((imag(alp1)>=0)-.5).*alp1;
%%
alp2  = sqrt(gam2.^2-el^2);
alp2  = sign((imag(alp2)>=0)-.5).*alp2;
%%
rts   = {alp1,alp2};

if do_test
   format long;
   for j=1:2
      h     = hh(j);
      E     = youngs(j);
      disp('test dispersion relation for T, h, E:');
      disp({T,h,E});
      %%
      H_dim    = H*L;%%reset in ND2 if too deep
      gam      = Rts{j}/L%%dimensional roots [m^{-1}]
      gam0     = [gam(1),2*pi/GEN_get_ice_wavelength(h,T,H_dim,E)]
      wavlen   = 2*pi./gam0
      %%
      om    = 2*pi/T;
      pram  = NDphyspram(0);
      g     = pram(2);
      rhow  = pram(3);
      rhoi  = pram(4);
      nu    = pram(5);
      %%
      D     = E*h^3/12/(1-nu^2);
      L5    = D/rhow/om^2;
      L1    = L5^.2;
      sig   = rhoi*h/rhow;
      alp   = om^2/g;
      c1    = 1/alp-sig;
      %%
      fdr   = disprel(gam*L1,1,c1/L1,H_dim/L1)
      if 1%%do a plot of real disp rel
         if 0
            kk = linspace(0,2*gam(1),51);
            ff = disprel(kk,L5,c1,H_dim);
         else
            kk = linspace(0,2*gam(1),51);
            ff = disprel(kk*L1,1,c1/L1,H_dim/L1);
            %%
            %gam0(2)*L1,1,c1/L1,H_dim/L1,L1
            f0 = disprel(gam0(2)*L1,1,c1/L1,H_dim/L1)
            f1 = disprel(gam0(1)*L1,1,c1/L1,H_dim/L1)
         end
         plot(kk,ff,'b',...
              kk,0*kk,'k',...
              gam0(2),0,'.r',...
              gam0(1),0,'or');
         ylim(5*[-1 1]);
         xlim(2*gam0(2)*[0 1]);
         pause;
      end
   end
end

function [f,df]   = disprel(gam,c5,c1,H)

if 1
   kt = gam.*tanh(gam*H);
   f  = (c5*gam.^4+c1).*kt-1;
   df = (5*c5*gam.^4+c1).*tanh(gam*H)+(c5*gam.^4+c1).*(gam*H).*sech(gam*H).^2;
else
   kt = gam.*sinh(gam*H);
   c  = cosh(gam*H);
   f  = (c5*gam.^4+c1).*kt-c;
   df = (5*c5*gam.^4+c1-H).*sinh(gam*H)+(c5*gam.^4+c1).*(gam*H).*c;
end
