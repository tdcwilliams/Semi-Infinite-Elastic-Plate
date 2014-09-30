function [x,y,z]=GEN_smooth2d(x,y,z);

%% smooth in x direction;
z  = .5*( z(3:end,:)+z(1:end-2,:) );
x  = x(2:end-1);


%% smooth in y direction;
z  = .5*( z(:,3:end)+z(:,1:end-2) );
y  = y(2:end-1);
