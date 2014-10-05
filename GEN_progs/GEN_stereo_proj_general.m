function [X,Y] = GEN_stereo_proj_general(lon,lat,lonc,latc)

DO_TEST  = 0;
if nargin==0
   lonc  = 45;
   latc  = 45;
   lat   = [-90:15:90]';
   lon   = lonc+0*lat;
   if 1
      lon   = [lon;lonc-30+0*lat;lonc+30+0*lat;lonc-60+0*lat;lonc+60+0*lat];
      lat   = [lat;lat;lat;lat;lat];
   end
   DO_TEST  = 1;
end

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

%% get orthogonal basis for the plane (xc,yc,zc)*(X,Y,Z)'=0;
[v1,v2]  = GEN_plane_basis([xc,yc,zc]');

%% find where the line from (x,y,z) to (-xc,-yc,-zc)
%%  intersects the plane (xc,yc,zc)*(X,Y,Z)'=0;

[xp,yp,zp]  = plane_intersect(x,y,z,xc,yc,zc);
[X,Y]       = convert_basis(v1,v2,xp,yp,zp);

if 0%% Keep projection as is;
   Mrot  = eye(2);
else
   %% projections of north and south poles 
   %% (rotate so that they both lie on the y axis);
   [xpNS,ypNS,zpNS]  = plane_intersect([0 0]',[0 0]',[1 -1]',xc,yc,zc);
   [X_NS,Y_NS]       = convert_basis(v1,v2,xpNS,ypNS,zpNS);
   if X_NS(1)~=Inf
      ang_NP   = pi/2;
      ang      = angle( X_NS(1)+1i*Y_NS(1) )-ang_NP;%% North pole should have ang==0
      Mrot     = [ [cos(-ang);sin(-ang)],[-sin(-ang);cos(-ang)] ];%%rotate vectors by -ang
   else
      ang_SP   = -pi/2;
      ang      = angle( X_NS(1)+1i*Y_NS(1) )-ang_SP;%% South pole should have ang==0
      Mrot     = [ [cos(-ang);sin(-ang)],[-sin(-ang);cos(-ang)] ];%%rotate vectors by -ang
   end
end

%% Rotate and convert to km;
XY = GEN_radius_earth*[Mrot(1,1)*X+Mrot(1,2)*Y,Mrot(2,1)*X+Mrot(2,2)*Y];
X  = XY(:,1);
Y  = XY(:,2);

if DO_TEST
   plot(X,Y,'o');
end

function [xp,yp,zp]  = plane_intersect(x,y,z,xc,yc,zc)

d1 = x+xc;
d2 = y+yc;
d3 = z+zc;
%%
t  = -(x*xc+y*yc+z*zc)./(xc*d1+yc*d2+zc*d3);
xp = x+d1.*t;
yp = y+d2.*t;
zp = z+d3.*t;

function [X,Y] = convert_basis(v1,v2,xp,yp,zp);

X  = v1(1)*xp+v1(2)*yp+v1(3)*zp;
Y  = v2(1)*xp+v2(2)*yp+v2(3)*zp;
