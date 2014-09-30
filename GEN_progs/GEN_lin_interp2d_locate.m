function [jx,jy,tx,ty] =...
            GEN_lin_interp2d_locate(x,y,xref,yref);

jx = 0*x;
tx = jx;
for j = 1:length(x)
   jx(j) = max(find(x(j)>=xref));
   x0    = xref(jx(j));
   if x~=x0
      x1    = xref(jx(j)+1);
      tx(j) = (x-x0)/(x1-x0);
   end
end

jy = 0*y;
ty = jy;
for j = 1:length(y)
%  y(j),yref,exp(y(j)),exp(yref)
%  find(y(j)>=yref)
   jy(j) = max(find(y(j)>=yref));
   y0    = yref(jy(j));
   if y~=y0
      y1    = yref(jy(j)+1);
      ty(j) = (y-y0)/(y1-y0);
   end
end
