function GEN_plot_stepfxn(x,y,colin)

nx = length(x);
X  = zeros(2*nx-2,1);
Y  = X;
%%
X(1:2:end)  = x(1:end-1); 
X(2:2:end)  = x(2:end);
Y(1:2:end)  = y(1:end-1); 
Y(2:2:end)  = y(1:end-1);
%%
if exist('colin')
   plot(X,Y,colin);
else
   plot(X,Y);
end
