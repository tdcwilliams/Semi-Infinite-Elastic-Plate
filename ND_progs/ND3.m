function y=ND3(Z,hh,N)
%% CALL: y=ND3(Z,[h1 h2]/h0,N);
%% Z={h,T,H_ice,th}, N=no of imaginary roots req'd;
%% output={K,{gam1,gam2},H,th,[lam,-mu],char,L}.

h=Z{1}; T=Z{2}; H_dim=Z{3}; thdeg=Z{4};
pram=NDphyspram(0);
H_infdep=5;%% H>H_infdep is the infinite depth criterion

E=pram(1);%% Pa
g=pram(2);%% m/s^2
rho=pram(3);%% kg/m^3
rho_ice=pram(4);%% kg/m^3
nu=pram(5);

om=2*pi/T;
D=E*h^3/12/(1-nu^2);
L_ice=(D/rho/g)^.25;
L5=D/rho/om^2; L=L5^.2;

lam=g/L/om^2; mu=rho_ice*h/rho/L; del=lam-mu;
H=H_dim/L; H=min(H,H_infdep);
K=RTS_ice_roots(del,H,N);

for j=1:2
	hr=hh(j);
    if hr==1
	Gam{j}=K;
    elseif hr~=0
	L2=hr^.6; H2=H/L2; del2=(lam-hr*mu)/L2;
	Gam{j}=RTS_ice_roots(del2,H2,N)/L2;
    else
	Gam{j}=RTS_wtr_roots(lam,H,N);
    end
end

char=[L_ice sqrt(L_ice/g)];
y={K,Gam,H,thdeg,[lam,-mu],char,L};
