function GEN_cheby_plot(f0,a,b,varargin);

linestyle='-';
if nargin==4
	linestyle=varargin{:};
end
%%
N=length(f0);
cheby_zeros=cos( pi/2/N*(2*(1:N)'-1) );
%%
nn=200; X=-1+(0:nn)'*2/nn;
x=a+(b-a)*(X+1)/2;
%%
Tn0=GEN_cheby_polys(cheby_zeros,N-1).';%{Tn0,f0}
Tn=GEN_cheby_polys(X,N-1);
%%
cn=2*Tn0*f0/N;
f=Tn*cn; f=f-(1+0*f)*diag(cn(1,:)/2);
plot(x,abs(f),linestyle);
