function [ag,wavlen]=GEN_get_wtr_groupvel(T,H)
%% CALL: [ag,wavlen]=GEN_get_wtr_wavelength(T,varargin)
%% calc's group velocity (ag) & wavelength ('wavlen') corresponding to a given
%% wave period 'T', and water depth
%% (optional argument - if no
%% value is specified infinite water depth is assumed)

if nargin==1
   H  = Inf;
end
[ag,wavlen]=GEN_get_ice_groupvel(0,T,H);
