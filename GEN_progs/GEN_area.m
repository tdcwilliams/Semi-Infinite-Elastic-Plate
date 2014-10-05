function Area=GEN_area(xy)
%% CALL: Area=GEN_area(xy)
%% Gives the area of polygon with vertices
%% if Area<1 going clockwise
%% xy = [x_1,y_1;
%%       x_2,y_2;
%%       .....
%%       x_N,y_N];

if nargin==0
   if 1
      %% test on area of a triangle:
      x3    = [-58.8993  -42.7814   26.0770]';
      y3    = [-15.1438   68.0290   38.3394]';
      xy    = [x3,y3];
      rs    = GEN_perimeter(xy,0);%%running sum
      p     = GEN_perimeter(xy,1);%%total perimeter
      abc   = [rs(2:end)-rs(1:end-1);p-rs(end)]
      s     = p/2;
      A     = sqrt(s*prod(s-abc))
   else
      x  = [-3 -1 6 3 -4]';
      y  = [-2 4 1 10 9]';
      xy = [x,y];
      A  = 60
   end

end

N     = size(xy,1);
xy    = [xy;xy(1,:)];
Area  = 0;

for j=1:N
   D     = det(xy(j:j+1,:))
   Area  = Area+D/2;
end
