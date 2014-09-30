function [gam,alp,H,el,L]=NDmakeZ_wtr(Z,N)
%% non-dimensionalization program for water problems when no ice is involved.
%% CALL: [gam,alp,H,el,L]=NDmakeZ_wtr(Z,N)
%% INPUTS: N is no of imag roots to use;
%%   Z=vector T or [T theta_inc=angle of incidence]
%%    or cell {T, theta_inc, water depth};
%% OUTPUTS: gam is vector of wave numbers;
%%   alp is component of gam in x direction;
%%   H is nondimensional water depth;
%%   el is component of incident wave number in y direction;
%%   L=g/(2*pi/T)^2 is natural length that we nondim wrt.

theta=0;
H_dim=1000;

if length(Z)==1
  T=Z;
elseif length(Z)==2
  T=Z(1);
  theta=Z(2);
elseif length(Z)==3
  T=Z{1};
  theta=Z{2};
  H_dim=Z{3};
end

%% some parameters:
g=9.81;
om=2*pi/T;
L=g/om^2;
H=H_dim/L;
lam=1;

%% get the roots:
gam=RTS_wtr_roots(lam,H,N);
el=gam(1)*sin(pi*theta/180);
alp=sqrt(gam.^2-el^2);%%will be pos real or pos imag