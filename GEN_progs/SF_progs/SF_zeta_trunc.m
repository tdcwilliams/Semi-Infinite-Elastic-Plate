function y=SF_zeta_trunc(z,N);

nvec=(1:N);
DO_TRANS=(size(z)*[-1;1]>0);
if DO_TRANS
  z=z.';
end
%
[NVEC,Z]=meshgrid(nvec,z);
y=zeta(z)-sum( NVEC.^(-Z),2 );
%%
if DO_TRANS
  y=y.';
end