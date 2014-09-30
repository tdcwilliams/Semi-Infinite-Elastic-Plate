function y=GEN_interp_exp(tt,An)

DO_TEST=0;

if iscell(An)==0
  Nterms=(length(An)-1)/2;
else
  Nterms=An{1};
end
nvec=-Nterms:Nterms;
Exp=exp(i*pi*tt*nvec);

if iscell(An)==0
  y=Exp*An;
else
  y=Exp;
end

if DO_TEST
  plot(tt,[real(y),imag(y)]);
end