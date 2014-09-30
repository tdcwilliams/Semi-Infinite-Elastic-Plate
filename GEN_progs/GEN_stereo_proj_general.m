function [X,Y] = GEN_stereo_proj_general(lon,lat,lonc,latc)

th    = pi*lon/180;
phi   = pi*lat/180;
z     = sin(phi);
x     = cos(phi).*cos(th);
y     = cos(phi).*sin(th);
%%
th    = pi*lonc/180;
phi   = pi*latc/180;
zc    = sin(phi);
xc    = cos(phi)*cos(th);
yc    = cos(phi)*sin(th);
%%
d1 = x+xc;
d2 = y+yc;
d3 = z+zc;

%% find where the line from (x,y,z) to (-xc,-yc,-zc)
%%  intersects the plane (xc,yc,zc)*(x,y,z)'=0;
t  = -(x*xc+y*yc+z*zc)./(xc*d1+yc*d2+zc*d3)
xp = x+d1*t;
yp = y+d2*t;
zp = z+d3*t;

%%get orthogonal basis for the plane (xc,yc,zc)*(x,y,z)'=0;
[v1,v2]  = GEN_plane_basis([xc,yc,zc]');
X        = v1(1)*xp+v1(2)*yp+v1(3)*zp;
Y        = v2(1)*xp+v2(2)*yp+v2(3)*zp;
