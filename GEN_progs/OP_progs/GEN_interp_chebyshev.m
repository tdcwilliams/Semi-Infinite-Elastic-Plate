function f=GEN_interp_chebyshev(tt,An)

disp('GEN_interp_chebyshev.m is old file; change to OP_*');

NgC=length(An)-1;
f=0*tt;

C0=1+f;
f=f+An(1)*C0;

if NgC==0
  return;
else
  C1=tt;
  f=f+An(2)*C1;
end

for its=2:NgC
  Cn=2*tt.*C1-C0;
  f=f+An(its+1)*Cn;
  C0=C1;
  C1=Cn;
end
