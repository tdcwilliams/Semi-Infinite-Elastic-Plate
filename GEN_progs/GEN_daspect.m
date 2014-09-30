function GEN_daspect(H,beta)

if nargin==1
   beta  = H;
   H     = gca;
end

X  = get(H,'XLim');
Y  = get(H,'YLim');
%%
dx = X(2)-X(1);
dy = Y(2)-Y(1);
%%
set(H,'DataAspectRatio',[beta/dy,1/dx,1]);
