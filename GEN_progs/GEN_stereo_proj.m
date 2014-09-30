function [x,y] = GEN_stereo_proj(lon,lat,chi);
%% CALL: [x,y] = GEN_stereo_proj(lon,lat,chi);
%% *chi is an angle that rotates the projection for better viewing;
%% *projection is from the south pole;

th    = pi*lat/180;
phi   = pi*lon/180;
zz    = sin(th);
rr0   = cos(th);
%%
r  = rr0./(1+zz);
x0 = r.*cos(phi);
y0 = r.*sin(phi);

%% now rotate by chi to make viewing more optimal;
%% chi=-pi/2 best for Fram Strait;
if ~exist('chi')
   chi   = -pi/2;
end
x     = x0*cos(chi)-y0*sin(chi);
y     = x0*sin(chi)+y0*cos(chi);
