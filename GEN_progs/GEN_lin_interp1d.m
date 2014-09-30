function y = GEN_lin_interp1d(x,xref,yref);

nr    = length(xref);
y     = 0*x;

%% INTERPOLATE TO GET UPPER VALUES OF x:
ju    = find(x>max(xref));
xu    = x(ju);
x0    = xref(end-1);
x1    = xref(end);
y0    = yref(end-1);
y1    = yref(end);
tx    = (xu-x0)/(x1-x0);
y(ju) = y0+tx*(y1-y0);

%% INTERPOLATE TO GET LOWER VALUES OF x:
jl = find(x<min(xref));
xl    = x(jl);
x0    = xref(1);
x1    = xref(2);
y0    = yref(1);
y1    = yref(2);
tx    = (xl-x0)/(x1-x0);
y(jl) = y0+tx*(y1-y0);

%% 'NORMAL' VALUES OF x:
jn    = find(x<=max(xref)&x>=min(xref));
xn    = x(jn);
nx    = length(x(jn));


if nx<=nr%% if x is shorter than xref,
         %% this way is quicker;
   for j = 1:nx
      J        = max(find(xn(j)>=xref));
      x0       = xref(J);
      y0       = yref(J);
      y(jn(j)) = y0;

      if xn(j)~=x0
         x1       = xref(J+1);
         y1       = yref(J+1);
         tx       = (xn(j)-x0)/(x1-x0);
         y(jn(j)) = y0+tx*(y1-y0);
      end
   end
else%% if xref is shorter than x,
    %% this way is quicker;
   J        = find(xn==xref(end));
   y(jn(J)) = yref(end);
   for j=1:nr-1
      x0       = xref(j);
      x1       = xref(j+1);
      y0       = yref(j);
      y1       = yref(j+1);
      J        = find(xn>=x0 & xn<x1);
      tx       = (xn(J)-x0)/(x1-x0);
      y(jn(J)) = y0+tx*(y1-y0);
   end
end
