function Bn=GEN_Bernoulli_nos(N)

do_test=0;

Bn=zeros(N+1,1);
Bn(1)=1;
if N>=1
  n=1;
  Bn(n+1)=-1/2;
end

for n=2:2:N
  r=n/2;
  mm=2*r+1;
  Bn(n+1)=-( Bn(1)+mm*Bn(2) )/(n+1);
  alp=1/(n+1);
  for k=1:r-1
    bet=(r+1-k)*(2*r+3-2*k)/k/(2*k-1);
    alp=bet*alp;
    Bn(n+1)=Bn(n+1)-alp*Bn(2*k+1);
  end
end

if do_test
  N2=22;
  jt=1:N+1;
  jt2=1:N2+1;
  Bn2=[1,-1/2,1/6,0,-1/30,0,1/42,0,-1/30,0,5/66,0,...
        -691/2730,0,7/6,0,-3617/510,0,43867/798,...
          0,-174611/330,0,854513/138]';
  if N<=N2
    jt2=jt;
  else
    jt=jt2;
  end
  tst=[Bn(jt),Bn2(jt2)]
end