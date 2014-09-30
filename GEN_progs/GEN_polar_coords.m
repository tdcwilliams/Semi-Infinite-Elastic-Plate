function [r,theta]=GEN_polar_coords(x,y)

z=x+i*y;
r=abs(z);
theta=angle(z);