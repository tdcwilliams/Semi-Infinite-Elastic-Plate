function T=GEN_get_wtr_period(wavlen,H)
%% CALL: T=GEN_get_wtr_period(wavlen,varargin)
%% calc's period corresponding to a given open water wavelength ('wavlen'),
%% given water depth (optional argument - if no
%% value is specified infinite water depth is assumed)

g  = NDphyspram(2);%% m/s^2
k  = 2*pi./wavlen;

if nargin==1%%infinite depth
  kt  = k;
else
  kt  = k.*tanh(k*H);
end

om = sqrt(g*kt);
T  = 2*pi./om;
