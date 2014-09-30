function [lon,lat] = GEN_stereo_invproj(x0,y0,chi);
%% CALL: [lon,lat] = GEN_stereo_invproj(x,y,chi);
%% *chi is angle that was included in the forward projection
%%  to make viewing more optimal;
%% *projection is from the south pole;

%% 1st rotate by -chi 
%% chi=-pi/2 best for Fram Strait;
if ~exist('chi')
   chi   = -pi/2;
end
if ~exist('x0')
   lon0     = [-10 0 30 0]';
   lat0     = [79 90 -30 -15]';
   ll0      = [lon0,lat0]
   [x0,y0]  = GEN_stereo_proj(lon0,lat0,chi);
end
x     = x0*cos(chi)+y0*sin(chi);
y     = -x0*sin(chi)+y0*cos(chi);
w     = x+1i*y;
r     = abs(w);
phi   = angle(w);
lon   = phi*180/pi;
%%
del   = atan(1./r);
th    = 2*del+pi/2;

%th       = 0*r;
%jj       = find(r<=1);
%th(jj)   = acos(r(jj));
%%%
%jj       = find(r>1);
%th(jj)   = -acos(1./r(jj));
lat      = 180/pi*th;
jj       = find(lat>90);
lat(jj)  = lat(jj)-180;
