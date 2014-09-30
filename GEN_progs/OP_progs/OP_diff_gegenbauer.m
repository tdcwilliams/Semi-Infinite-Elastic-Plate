function [Md,alp_out]=OP_diff_gegenbauer(alp_in,Ngg,m,Lc)
%% CALL: [Md,alp_out]=OP_diff_gegenbauer(alp_in,Ngg,m,Lc)
%% INPUTS:
%% alp_in & Ngg are the gegenbauer index and polynomial degree
%%  of the series to be differentiated -
%%  ie want to diff \sum_{n=0}^Ngg A_nC_n^{alp_in};
%% m is the no of times to differentiate;
%% LC is the half-length of the interval
%%  that the approximation is over;
%% OUTPUTS:
%% alp_out is gegenbauer index of new series;
%% Md*An gives the coefficients of the new series.
%%
%% NB can let Ngg=cell{ An };
%% then output Md changes to Md*An
%% ie new coefficients are outputted.
%%
%% C'_n^(alpha) = 2\alpha*C_{n-1}^(alpha+1)

want_coeffs=iscell(Ngg);
if want_coeffs
  An=Ngg{1};
  Ngg=length(An)-1;
end

if m==0
  if want_coeffs
    Md=An;
  else
    Md=eye(Ngg+1);
  end
  alp_out=alp_in;
  return;
end

nrows=Ngg+1-m;
zz=zeros(nrows,m);
alp_out=alp_in+m;

if alp_in~=0%%C_n^(alpha)->2\alpha*C_{n-1}^(alpha+1)
  Md=eye(nrows);
  fac=2^m*gamma(alp_in+m)/gamma(alp_in);
  Md=[zz,fac*Md];
else%%T_n'=nU_{n-1}
  Md=diag(m:Ngg);
  fac=2^(m-1)*gamma(alp_in+m)/gamma(alp_in+1);
  Md=[zz,fac*Md];
end

if nargin==4
  Md=Md/Lc^m;
end

if want_coeffs
  Md=Md*An;
end