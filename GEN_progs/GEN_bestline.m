function [m,c]=GEN_bestline(x,y,czero)

if nargin==2 | czero==0
  mc  = polyfit(x,y,1);
  m   = mc(1);
  c   = mc(2);
else%% require c=0
  c   = 0;
  m   = (x'*y)/(x'*x);
end
