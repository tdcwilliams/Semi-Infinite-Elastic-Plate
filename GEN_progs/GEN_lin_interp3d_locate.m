function [jx,jy,jz,tx,ty,tz] =...
            GEN_lin_interp3d_locate(x,y,xref,yref);

jx = 0*x;
tx = jx;
for j = 1:length(x)
   jx(j) = max(find(x(j)>=xref));
   x0    = xref(jx(j));
   if x~=x0
      x1    = xref(jx(j)+1);
      tx(j) = (x(j)-x0)/(x1-x0);
   end
end

jy = 0*y;
ty = jy;
for j = 1:length(y)
   jy(j) = max(find(y(j)>=yref));
   y0    = yref(jy(j));
   if y~=y0
      y1    = yref(jy(j)+1);
      ty(j) = (y(j)-y0)/(y1-y0);
   end
end

jz = 0*z;
tz = jz;
for j = 1:length(z)
   jz(j) = max(find(z(j)>=zref));
   z0    = zref(jz(j));
   if z~=z0
      z1    = zref(jz(j)+1);
      tz(j) = (z(j)-z0)/(z1-z0);
   end
end
