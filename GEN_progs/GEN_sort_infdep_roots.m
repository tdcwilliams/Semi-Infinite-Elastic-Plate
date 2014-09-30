function gam=GEN_sort_infdep_roots(r);
%% CALL: gam=GEN_sort_infdep_roots(r);
%% r=roots([1 0 0 0 del -1]);
%% gam(1)>0 is real;
%% gam(2)=conj(gam(3)) is in 1st quadrant;
%% gam(4:5) are in lh half-plane.

gam=0*r;
gam(1)=r( find(r>0 & imag(r)==0) );
gam(2:3)=r( find(r>0 & imag(r)~=0) );
gam(4:5)=r (find(r<0) );
