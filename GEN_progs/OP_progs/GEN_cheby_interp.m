function [f,x]=GEN_cheby_interp(f0,a,b);

N=length(f0);
cheby_zeros=cos( pi/2/N*(2*(1:N)'-1) );
%%
nn=1000; X=-1+(0:nn)'*2/nn;
x=a+(b-a)*(X+1)/2;
%%
Tn0=GEN_cheby_polys(cheby_zeros,N-1).';
Tn=GEN_cheby_polys(X,N-1);%{Tn,Tn0,f0}
%%
cn=2*Tn0*f0/N;
f=Tn*cn; f=f-(1+0*f)*diag(cn(1,:)/2);