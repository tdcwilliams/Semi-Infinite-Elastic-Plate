function M=GEN_swap_chebyshev(Np,dir)

if dir==0% Matrix changes Tn coeff's to Un coeff's:
  M=.5*eye(Np+1);
  jj=1:Np-1;
  M(jj+2,jj)=M(jj+2,jj)-M(jj,jj);
  M(1)=1;
else% Matrix turns Un coeff's to Tn coeff's:
  D=2*eye(Np+1);
  D(1)=1;
  n0=3;
  N=Np-1;
  M=D;
  while n0<=Np+1
    jj=n0:Np+1;
    M(jj,1:N)=M(jj,1:N)+D(1:N,1:N);
    n0=n0+2;
    N=N-2;
  end
end
