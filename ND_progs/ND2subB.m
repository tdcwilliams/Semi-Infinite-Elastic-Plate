function [Rts,HHsig,th,del0,L]=ND2subB(Z,hh,Nmax);
%% CALL: y=ND2(Z,hr,N);
%% Z={h,T,H_dim,th}, hr=h1/h0, N=no of imaginary roots req'd;
%% output={K,gam,H,th,[lam,-mu],char,L}.

T=Z{1}; thdeg=Z{2}; th=thdeg*pi/180; H_dim=Z{3};
pram=NDphyspram(0);
H_infdep=5;%% H>H_infdep is the infinite depth criterion

E=pram(1);%% Pa
g=pram(2);%% m/s^2
rho=pram(3);%% kg/m^3
rho_ice=pram(4);%% kg/m^3
nu=pram(5);

hmax=max(hh); hr=hh/hmax; Dr=hr.^3;

om=2*pi/T;
D=E*hmax^3/12/(1-nu^2);
L_ice=(D/rho/g)^.25;
L=(D/rho/om^2)^.2;

lam=g/L/om^2; sigsig=rho_ice*hh/rho/L; H=H_dim/L;
sig=abs(sigsig*[1;-1]); HH=H-sigsig;
dH=max(HH)-H_infdep; HH=HH-dH*(dH>0);% [HH(2)-HH(1), sig]
HHsig=[HH,sig]; mu=max(sigsig); del0=[lam,-mu];

Nsig=round(sig/max(HH)*Nmax); %N1=Nmax-Nsig-3;
%Nmax=NN(1); N_sig=NN(2);
if hh(1)==hh(2)
	gam=RTS_ice_roots(lam-mu,HH(1),Nmax);
	Rts={gam,gam,[]};
else
	jmax=find(hh==hmax); jmin=find(hh~=hmax);
	Rts{jmax}=RTS_ice_roots(lam-mu,HH(jmax),Nmax);
    if Dr(jmin)==0
	Rts{jmin}=RTS_wtr_roots(lam,HH(jmin),Nmax+Nsig+1);
	Rts{3}=RTS_wtr_roots(lam,sig,Nsig);
    else
	L_j=Dr(jmin)^.2; del_j=(lam-sigsig(jmin))/L_j;
	Rts{jmin}=RTS_ice_roots(del_j,HH(jmin)/L_j,Nmax+Nsig+3)/L_j;
	Rts{3}=RTS_ice_roots(del_j,sig/L_j,N_sig)/L_j;
    end
end
