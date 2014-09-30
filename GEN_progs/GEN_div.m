function [y,rem]=GEN_div(x,base);
y0=x/base; y=round(y0-.5);
y=y+(y==(y0-1)); rem=x-y*base;
