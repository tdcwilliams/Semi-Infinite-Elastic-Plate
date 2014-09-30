%% GEN_Principal_Axes.m
%% Author: Timothy Williams
%% Date:   20130619, 08:39:24 CEST

function [xy_rel,xy_av,U] = GEN_Principal_Axes(xy)

DO_TEST  = 0;
if nargin==0
   DO_TEST  = 1;

   %%Unrotated ellipse for testing;
   a     = 4;%%make a>b so R(:,1) (rotation of x axis) is major axis
   b     = 2;
   th    = linspace(0,2*pi,1000)';
   x0    = a*cos(th);
   y0    = b*sin(th);

   %%Rotate ellipse:
   rot   = pi/6;
   R     = [cos(rot) -sin(rot);sin(rot) cos(rot)];
   xy    = [x0 y0]*R';

   %%Translate ellipse:
   x1 = .2;
   y1 = .4;
   xy(:,1)  = xy(:,1)+x1;
   xy(:,2)  = xy(:,2)+y1;

   %%Correct principal axes;
   r1 = R(:,1)
   r2 = R(:,2)
end




%% subtract mean & find principal values:
%% xy = [x y]: N x 2 matrix
xy_av = mean(xy);
for j=1:2
   xy(:,j)  = xy(:,j)-xy_av(j);
end
A           = xy'*xy;%%(2 x N).(N x 2) = (2 x 2) 
[U,Lam,V]   = svd(A);
%%
xy_rel      = xy*U;%%(U'*xy')'


if DO_TEST==1
   x  = xy(:,1);
   y  = xy(:,2);
   plot(x,y)
   %%
   hold on;
   daspect([1 1 1]);
   X  = 1.1*max(x);
   Y  = 1.1*max(y);
   axis([ X*[-1 1],Y*[-1 1] ]);
   %%
   plot(X*[-1 1],[0 0],'k')
   plot([0 0],Y*[-1 1],'k')
   %%
   rr    = {r1,r2};
   colin = {'--r','--g'};
   for j=1:2
      vv = rr{j};
      xx = 1.1*max(a,b)*vv(1)*[-1 1];
      yy = 1.1*max(a,b)*vv(2)*[-1 1];
      plot(xx,yy,colin{j});
   end

   hold off;

   figure(2);
   x  = xy_rel(:,1);
   y  = xy_rel(:,2);
   plot(x,y)
   daspect([1 1 1]);

end
