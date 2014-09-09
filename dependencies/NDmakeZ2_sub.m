function [Rts,rts,HHsig,el,del0,L,rts_approx]=...
            NDmakeZ2_sub(Z,hh,NNN,EE,rho)
%% CALL: [Rts,rts,HHsig,el,del0,L,rts_approx]
%% =NDmakeZ2_sub(Z,hh,NNN,youngs)
%% Z=T or [T th] or {T th H_dim}

do_test  = 0;
if nargin==0%%do a test
   do_test  = 1;
   hh       = [1 2];
   %youngs   = [5e9 6e9];
   %youngs   = [5e9 5e9];
   youngs   = NDphyspram(1)*[1 1];
   rho      = 1025;
   NNN      = 10;
   %%
   T     = 5;
   theta = 0;
   H_dim = 100;
   Z     = {T,theta,H_dim};
end

%%sometimes need to solve a 3rd dispersion relation, but sometimes don't...
if iscell(NNN)
   WANT3 = (length(NNN)==4);
else
   WANT3 = (length(NNN)==3);
   if length(NNN)==1
      NNN   = [NNN,NNN];
   end
end

%% Z->{T,th,H_dim};
if length(Z)==1
   T        = Z;
   Z_       = NDphyspram([]);
   Z        = {T};
   Z(2:3)   = Z_([4 3]);
elseif length(Z)==2
   T     = Z(1);
   theta = Z(2);
   Z_    = NDphyspram([]);
   Z     = {T,theta,Z_{3}};
%  elseif length(Z)==3
%    Z(2:3)=Z([3 2]);
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

%% get impt parameters:
[Rts,HHsig,th,del0,L]   = ND2sub(Z,hh,NNN,EE,rho);

%[Rts{1}(1),Rts{2}(1)],pause

%% sort out roots:
gam0  = Rts{1}(1);
el    = gam0*sin(th);
for j=1:2+WANT3
   gam      = Rts{j};
   alp      = sqrt(gam.^2-el^2);
   rts{j}   = sign((imag(alp)>=0)-.5).*alp;
   %rts_approx{j}=( 1:NNN(j) )'*pi*i/HHsig(j);
end

%  if HHsig(3)==0
%    rts_approx{3}=[];
%  end

if do_test
   format long;
   for j=1:2
      h     = hh(j);
      E     = youngs(j);
      disp('test dispersion relation for T, h, E:');
      disp({T,h,E});
      %%
      gam   = Rts{j}/L%%dimensional roots [m^{-1}]
      gam0  = [gam(1),2*pi/GEN_get_ice_wavelength(h,T,H_dim,E)]
      %%
      om    = 2*pi/T;
      pram  = NDphyspram(0);
      rhow  = pram(2);
      rhoi  = pram(3);
      g     = pram(4);
      nu    = pram(5);
      %%
      D  = E*h^3/12/(1-nu^2);
      c5 = D/rhow/om^2;
      c1 = g/om^2-rhoi*h/rhow;
      kt = gam.*tanh(gam*H_dim);
      %%
      fdr   = (c5*gam.^4+c1).*kt%%should be 1
   end
end
