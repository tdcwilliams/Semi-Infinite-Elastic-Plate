function Y=GEN_setzero(x,Y,Xints)
%% CALL: GEN_setzero(x,y,xints)
%% if y is a function st y>=0
%% this sets the min value in the interval x0<x<x1
%% of y to 0.

N=size(Y,2);
for n=1:N
  y=Y(:,n);
  xints=Xints{n};%,n
  for j=1:size(xints,1)
    x0=xints(j,1);
    x1=xints(j,2);
    jj=find((x>x0) & (x<x1));
    yj=y(jj);
    jmin=find(yj==min(yj));%jmin
    yj(jmin(1))=0;
    y(jj)=yj;
  end
  Y(:,n)=y;
end