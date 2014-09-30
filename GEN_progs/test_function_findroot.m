function [y,dy]=test_function_findroot(x)
%  y=cos(x)-1;
a=1e-12;
na=2;
b=.3;
nb=3;
y=(x-a).^na.*sin(x).*(x-b).^nb;
dy=na*(x-a).^(na-1).*sin(x).*(x-b).^nb+...
    +(x-a).^na.*cos(x).*(x-b).^nb+...
    +nb*(x-a).^na.*sin(x).*(x-b).^(nb-1);