function y=testroots(period,h,H_dim)

addpath ../ND_progs

th=0;
Z=[period,th,h];
[del0,L]=NDbasics(Z,1);

del=sum(del0);
H=H_dim/L;

y=RTS_ice_roots(del,H,5)/L;