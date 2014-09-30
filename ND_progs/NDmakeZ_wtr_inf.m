function [gam,alp,el,L]=NDmakeZ_wtr_inf(Z,N)
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

th=0;

if length(Z)==1
  T=Z;
elseif length(Z)==2
  T=Z(1);
  th=pi*Z(2)/180;
end

%% some parameters:
g=9.81;
om=2*pi/T;
L=g/om^2;
gam=1;
el=sin(th);
alp=cos(th);