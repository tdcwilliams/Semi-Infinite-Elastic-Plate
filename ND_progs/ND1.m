function y=ND1(Z,N)
%CALL: y=Nondim(Z,N);
%INPUT:
%Z={h,T,H_dim,th} or {K,aHfac,th,del0,char,L};
%N=no of roots to use in e/fxn exp of V_{z,zeta}
%OUTPUT:
%{K,aHfac,th,[lam,-mu],char,L} (del=lam-mu).

h=Z{1};
T=Z{2};
H_dim=Z{3};
theta=Z{4};
inf_crit=1;
H_inf=5;
pram=NDphyspram(0);

E=pram(1); %Pa
g=pram(2); %m/s^2
rho=pram(3); %kg/m^3
rho_ice=pram(4); %kg/m^3
nu=pram(5);

om=2*pi/T;
D=E*h^3/12/(1-nu^2);
L_ice=(D/rho/g)^.25;
L5=D/rho/om^2;
L=L5^.2;

lam=g/L/om^2;
mu=rho_ice*h/rho/L;
del=lam-mu;
H=H_dim/L;
if inf_crit==1 %APPLY INFINITE DEPTH CRITERION
  H=min(H_inf,H);
  aH(2)=H;
end
K=RTS_ice_roots(del,H,N);%%changed from original
char=[L_ice sqrt(L_ice/g)];
y={K,H,theta,[lam,-mu],char,L};