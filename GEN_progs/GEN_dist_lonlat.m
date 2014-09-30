function s  = GEN_dist_lonlat(ll1,ll2)
%% CALL: s  = GEN_dist_lonlat(ll1,ll2)
%% s is a vector with the distances (in km) between
%%  the points in the rows in ll1 & ll2;
%%       ll1=[lon_1, lat_1; lon_2,lat_2;...];
%% 
%% CALL: s  = GEN_dist_lonlat(ll1)
%% s is a vector containing the successive distances (in km)
%%  between the points in the rows of ll1;
%%
%% if ll1 is a transect, and you want its total length S,
%% CALL: S=sum( GEN_dist_lonlat(ll1) );

if nargin==0
   if 0
      ll1   = [0 79];
      ll2   = [1 79];
   else
      ll1   = [0 79;1 79];
   end
   [dx,radius_earth]=GEN_lon2distance(ll1(1,2))
end

s  = 0;
if exist('ll2')
   phi1     = pi/180*ll1(:,2);%%lat
   phi2     = pi/180*ll2(:,2);%%lat
   %%
   th1   = pi/180*ll1(:,1);
   th2   = pi/180*ll2(:,1);
else
   if size(ll1,1)==1%%only one point;
      'hi'
      return;
   end
   %%
   phi1  = pi/180*ll1(1:end-1,2);
   phi2  = pi/180*ll1(2:end,2);
   %%
   th1   = pi/180*ll1(1:end-1,1);
   th2   = pi/180*ll1(2:end,1);
end

dphi     = phi2-phi1;
phi_av   = (phi2+phi1)/2;
dth      = th2-th1;
%%
r  = 6371;%%radius of earth (km);
s  = r*sqrt(dphi.^2+cos(phi_av).^2.*dth.^2);
