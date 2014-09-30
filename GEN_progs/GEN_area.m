function Area=GEN_area(xy)
%% CALL: Area=GEN_area(xy)
%% Gives the area of polygon with vertices
%% if Area<1 going clockwise
%% xy = [x_1,y_1;
%%       x_2,y_2;
%%       .....
%%       x_N,y_N];

N     = size(xy,1);
xy    = [xy;xy(1,:)];
Area  = 0;

for j=1:N
  Area   = Area+abs(det(xy(j:j+1,:))/2);
end
