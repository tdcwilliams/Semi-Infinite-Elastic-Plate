function [w,del]=plot_roots(T,H)
%% call: [w,del]=plot_roots(T,H);
%% T is the period, H the nondimensional water depth;
%% del is the main parameter of the dispersion relation:
%% f=(k^4+del)*tanh(k*H)-1=0; w=k*H/pi (ie is the roots scaled).

[y del0]=prams(1,T); del=sum(del0);
w=RTS_ice_roots(del,H,2)*H/pi;
plot(w,'rx'), hold on, plot(-w,'rx');
X=max(real(w))*1.2*[-1 1]; Y=max(imag(w))*1.2*[-1 1];
plot(X,0*X,':'); plot(0*Y,Y,':'); hold off;
axis([X Y]);

function [tau_bar,del0,y,wavlen]=prams(h,T)
%% [tau_bar,del0,y,wavlen]=prams(h,T); T can be a vector,
%% tau-bar is the nondimensional period;
%% wavlen is the nondimensional, infinite depth wavelength
%% for an ice-coupled flexural-%%gravity wave; y={L_ice,L,sigma(h\,m)}
%% (L_ice is the characteristic length, L is the natural length, 
%% sigma is a surface density type parameter).
pram=physpram(0);
E=pram(1); %Pa
g=pram(2); %m/s^2
rho=pram(3); %kg/m^3
rho_ice=pram(4); %kg/m^3
nu=pram(5);

om=2*pi./T;
D=E*h^3/12/(1-nu^2);
L_ice=(D/rho/g)^.25; T_ice=sqrt(L_ice/g);
L5=D/rho./om.^2; L=L5.^.2;
lam=g./L./om.^2; mu=rho_ice*h/rho./L; sigma=mu.*L/L_ice;
tau_bar=T/T_ice; del0=[lam,-mu]; y={L_ice,L,sigma(1)};

if nargout==4
	del=lam-mu;
    for j=1:length(T)
	r=roots([1 0 0 0 del(j) -1]);
	gam0=r( find(imag(r)==0 & r>0) ); wavlen(j,1)=2*pi/gam0;
    end
end

function y=physpram(n);
%%%stores the default physical parameters to use
%% CALL: y=physpram(n)
%% INPUT:
%%    n=[]:->y={h,T,aH_dim,th,[L0,T0,m0]};
%%    n=0:->y=pram=[E g rho rho_ice nu];
%%    else, y=pram(n).
E=5e9;%Pa
g=9.81;%m/s^2
rho=1025;%kg/m^3
rho_ice=922.5;%kg/m^3
nu=.3;

if isempty(n)==1
	h=1; T=5; H_dim=70; th=0;
	D=E*h^3/12/(1-nu^2); 
	L0=(D/rho/g)^.25; T0=sqrt(L0/g); m0=rho_ice*h/rho/L0;
	y={h,T,H_dim,th,[L0,T0,m0]};
else
	y0=[E g rho rho_ice nu];
	if n==0
		y=y0;
	else
		y=y0(n);
	end
end
