function [jx,tx] =...
            GEN_lin_interp1d_locate(x,xref);

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
