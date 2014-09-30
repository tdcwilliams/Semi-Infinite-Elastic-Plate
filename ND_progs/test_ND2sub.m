Z     = {10,0,100};
E0    = 5.45e9;
hh    = [1 2];
NNN   = [10,10];
rho   = 1025;


EE                      = [E0,E0;.9*rho,.9*rho;.3,.3];
[Rts,HHsig,th,del0,L]   = ND2sub(Z,hh,NNN,EE,rho);
k11                     = Rts{1}/L;
k21                     = Rts{2}/L;

EE                      = [2*E0,2*E0;.9*rho,.9*rho;.3,.3];
[Rts,HHsig,th,del0,L]   = ND2sub(Z,hh,NNN,EE,rho);
k12                     = Rts{1}/L;
k22                     = Rts{2}/L;

EE                      = [1*E0,2*E0;.9*rho,.9*rho;.3,.3];
[Rts,HHsig,th,del0,L]   = ND2sub(Z,hh,NNN,EE,rho);
k13                     = Rts{1}/L;
k23                     = Rts{2}/L;

EE                      = [2*E0,1*E0;.9*rho,.9*rho;.3,.3];
[Rts,HHsig,th,del0,L]   = ND2sub(Z,hh,NNN,EE,rho);
k14                     = Rts{1}/L;
k24                     = Rts{2}/L;

[k11,k13]
[k12,k14]
[k21,k24]
[k22,k23]

EE                      = 1*E0;
[Rts,HHsig,th,del0,L]   = ND2sub(Z,hh,NNN,EE,rho);
k15                     = Rts{1}/L;
k25                     = Rts{2}/L;
[k11,k15]
[k21,k25]

EE                      = 2*E0;
[Rts,HHsig,th,del0,L]   = ND2sub(Z,hh,NNN,EE,rho);
k16                     = Rts{1}/L;
k26                     = Rts{2}/L;
[k12,k16]
[k22,k26]
