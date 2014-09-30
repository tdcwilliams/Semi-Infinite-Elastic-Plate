function y=GEN_split_ri(x)

sz=size(x);
if sz(2)==length(x);%%row vector
  y=[real(x);imag(x)];
else%%column vector
  y=[real(x),imag(x)];
end