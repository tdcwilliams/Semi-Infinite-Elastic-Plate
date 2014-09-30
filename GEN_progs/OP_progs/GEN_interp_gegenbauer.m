function f=GEN_interp_gegenbauer(tt,An,alpha)

disp('please use OP_interp_gegenbauer.m');

NgC=length(An)-1;
f=0*tt;

C0=f;
C1=1+f;
f=f+An(1)*C1;

for its=1:NgC
  Cn=2*(its+alpha-1)/its*tt.*C1-(its+2*alpha-2)/its*C0;
  f=f+An(its+1)*Cn;
  C0=C1;
  C1=Cn;
end
