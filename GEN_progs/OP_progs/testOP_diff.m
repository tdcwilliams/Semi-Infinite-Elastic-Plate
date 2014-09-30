function testOP_diff(alp_in,Ngg,m)

k=3; Lc=2;
np=100;
tt=-1+2*(0:np)'/np;
xx=Lc*tt;

if alp_in~=0%% test OP_diff_gegenbauer.m
  [tgg,wgg]=OP_numint_jacobi(alp_in-.5,alp_in-.5,Ngg);
  ipC=OP_inprod_gegenbauer(tgg,wgg,alp_in);
  ygg=exp(i*k*Lc*tgg);
  an=ipC*ygg;
  %%
  dy=(i*k)^m*exp(i*k*xx);
  Md=OP_diff_gegenbauer(alp_in,Ngg,m,Lc);
  An=Md*an;
  dyap=OP_interp_gegenbauer(tt,alp_in+m,An);
  plot(xx,[real(dy),imag(dy)]); hold on;
  plot(xx,[real(dyap),imag(dyap)],'--r'); hold off;
 return
elseif 0%% test OP_diff_gegenbauer.m when alp_in==0
  [tgg,wgg]=OP_numint_chebyshev(Ngg);
  ipC=OP_inprod_chebyshev(tgg,wgg,Ngg);
  ygg=exp(i*k*Lc*tgg);
  an=ipC*ygg;
  %%
  dy=(i*k)^m*exp(i*k*xx);
  Md=OP_diff_gegenbauer(alp_in,Ngg,m,Lc);
  An=Md*an;
  if m==0
    dyap=OP_interp_chebyshev(tt,An);
  else
    dyap=OP_interp_gegenbauer(tt,alp_in+m,An);
  end
  plot(xx,[real(dy),imag(dy)]); hold on;
  plot(xx,[real(dyap),imag(dyap)],'--r'); hold off;
  return
else
  [tgg,wgg]=OP_numint_chebyshev(Ngg);
  ipC=OP_inprod_chebyshev(tgg,wgg,Ngg);
  ygg=exp(i*k*Lc*tgg);
  an=ipC*ygg;
  %%
  dy=(i*k)^m*exp(i*k*xx);
  Md=eye(Ngg+1);
  Nd=Ngg;
  for j=1:m
    Md0=OP_diff_chebyshev(Nd,Lc);
    Nd=Nd-1;
    Md=Md0*Md;
  end
  An=Md*an;
  dyap=OP_interp_chebyshev(tt,An);
  plot(xx,[real(dy),imag(dy)]); hold on;
  plot(xx,[real(dyap),imag(dyap)],'--r'); hold off;
  return
end