function [dAn,alp_out,bet_out]...
		=OP_diff_jacobi(alp_in,bet_in,An,m,Lc)

want_coeffs=~iscell(An);
if ~want_coeffs%%An cell=>want matrix
  Ngg=An{1};
  An=ones(Ngg+1,1);
else%%get coeffs
  Ngg=length(An)-1;
end

alp_out=alp_in+m;
bet_out=bet_in+m;
Nout=Ngg-m;

%% DIFFERENTIATE m TIMES, USING
%% d/dx.P^(alf,bet)_n=(n+alf+bet+1)/2*P^(alf+1,bet+1)_{n-1}
%% (MATHWORLD):
n_need=(m:Ngg)';
dAn=An(n_need+1);
for j=1:m
  fac=(n_need+alp_in+bet_in+j)/2;
  dAn=dAn.*fac;
end

%% ALLOW FOR INTERVALS OTHER THAN [-1,1]:
if nargin==5
  dAn=dAn/Lc^m;
end

%% IF HAVEN'T PROVIDED A SPECIFIC VECTOR OF COEFF'S,
%% MAKE A MATRIX THAT'LL DIFFERENTIATE AN ARBITRARY ONE:
if ~want_coeffs
  zz=zeros(Nout+1,Ngg-Nout);
  dAn=[zz,diag(dAn)];
end

if 0%%TEST ROUTINE ON e^{ikx} (TESTED & OK, 16/7/07):
  k=3;
  if want_coeffs
    zz=zeros(Nout+1,Ngg-Nout);
    dAn=[zz,diag(dAn)];
  end
  [tgl,wgl]=OP_numint_jacobi(alp_in,bet_in,Ngg);
  ipJ=OP_inprod_jacobi(tgl,wgl,alp_in,bet_in,Ngg);
  fj=real(exp(i*k*tgl));
  fn=ipJ*fj;
  %%
  np=100;
  x=-1+2*(0:np)'/np;
  df=real( (i*k)^m*exp(i*k*x) );
  plot(x,df), hold on;
  %%
  dfn=dAn*fn;
  df2=OP_interp_jacobi(x,alp_out,bet_out,dfn);
  plot(x,df2,'--r'), hold off;
end
