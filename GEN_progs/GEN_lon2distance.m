function [dx,radius_earth]=GEN_lon2distance(lat)
%% CALL: [dx,radius_earth]=GEN_lon2distance(lat)
%% dx is distance in km of 1 degree of longitude at latitude 'lat';
%% radius_earth is in km and is the mean radius;

radius_earth   = 6371;%%km
dx             = (pi/180)*radius_earth*cos(pi*lat/180);
