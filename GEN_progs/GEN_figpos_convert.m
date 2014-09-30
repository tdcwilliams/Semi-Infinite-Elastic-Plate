function [X,Y] = GEN_figpos_convert(x,y)
%% CALL: [X,Y] = GEN_figpos_convert(x,y)
%% GEN_figpos_convert.m
%% Author: Timothy Williams
%% Date:   20130626, 09:54:19 CEST
%% converts (x,y) to values between 0,1 - relative to xlim,ylim values

xl = get(gca,'xlim');
yl = get(gca,'ylim');
x_nan = ((x<xl(1))&(x>xl(2)));
y_nan = ((y<yl(1))&(y>yl(2)));

if x_nan|y_nan
   X  = NaN;
   Y  = NaN;
   return;
end

%%coord's relative to gca
Xo = (x-xl(1))/(xl(2)-xl(1));
Yo = (y-yl(1))/(yl(2)-yl(1));

%%coord's relative to gca
ap = get(gca,'position');
X  = Xo*ap(3)+ap(1);
Y  = Yo*ap(4)+ap(2);

return

%% example plot
clf
plot(rand(5,2)*5)
%% get info specific to the axes you plan to plot into
set(gcf,'Units','normalized')
set(gca,'Units','normalized')
ax = axis;
ap = get(gca,'Position')
disp( 'annotation from (1,2) to (3,4)')
xo = [1,3];
yo = [2,4];

hold on;
plot(xo(1),yo(1),'xk');
plot(xo(2),yo(2),'xk');
hold off;

Xo = (xo-ax(1))/(ax(2)-ax(1))
Yo = (yo-ax(3))/(ax(4)-ax(3))


%xp = Xo*(ap(3)-ap(1))+ap(1)
%yp = Yo*(ap(4)-ap(2))+ap(2)
xp = Xo*ap(3)+ap(1)
yp = Yo*ap(4)+ap(2)
ah=annotation('arrow',xp,yp,'Color','r');
