function GEN_plot_matrix(x,y,M,Mlims)

if nargin==3
   m0 = min(min(M));
   m1 = max(max(M));
else
   m0 = Mlims(1);
   m1 = Mlims(2);
end
%%
nx = length(x);
ny = length(y);
lims  = [ min(x),max(x),min(y),max(y) ];
%%
if length(x)==1&length(y)==1
   dx = 1;
   dy = 1;
elseif length(x)==1
   dy = y(2)-y(1);
   dx = dy;
elseif length(y)==1
   dx = x(2)-x(1);
   dy = dx;
else
   dx = x(2)-x(1);
   dy = y(2)-y(1);
end
dx0   = dx/2;
dy0   = dy/2;
%%
for i=1:ny
   for j=1:nx
      %% vertices of grid cell;
      X     = [x(j)-dx0,x(j)-dx0,x(j)+dx0,x(j)+dx0];
      Y     = [y(i)-dy0,y(i)+dy0,y(i)+dy0,y(i)-dy0];
      C     = (M(i,j)-m0)/(m1-m0);%% colour;
      pp    = patch(X,Y,C); hold on;
      if 0%max(nx,ny)>60
         set(pp,'LineStyle','none');
      end
      axis(lims);
   end
end
