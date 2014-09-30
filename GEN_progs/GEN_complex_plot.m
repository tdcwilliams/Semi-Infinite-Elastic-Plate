function GEN_complex_plot(zvec,char);

if nargin<2
  char='-k';
end
x=real(zvec);
y=imag(zvec);
plot(x,y,char);

