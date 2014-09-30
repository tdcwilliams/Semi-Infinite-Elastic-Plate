function Fn=GEN_fibonacci_nos(N)

do_test=0;

Fn=zeros(N+1,1);
if N>=1
  Fn(2)=1;
end
for n=2:N
  Fn(n+1)=Fn(n)+Fn(n-1);
end

if do_test
  N2=8;
  jt=1:N+1;
  jt2=1:N2+1;
  Fn2=[0 1 1 2 3 5 8 13 21]';
  if N<=N2
    jt2=jt;
  else
    jt=jt2;
  end
  tst=[Fn(jt),Fn2(jt2)]
end