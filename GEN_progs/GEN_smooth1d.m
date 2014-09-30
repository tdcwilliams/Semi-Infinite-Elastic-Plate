function [x,y]=GEN_smooth1d(x,y);

y  = .5*( y(3:end,:)+y(1:end-2,:) );
x  = x(2:end-1,:);
