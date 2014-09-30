function [gam1,H,gam0vec,del0,L,char]=ND1_outer3d(Z,N,varargin)
%% CALL: [gam1,H,gam0vec,del0,L,char]=ND1_outer3d(Z,N)
%% 
%% Z=T, [T th], {T,th,H_dim} or {h,T,H_dim,th}
%% N is the number of imaginary roots wanted.
%% 
%% OUTPUT: 
%% H is the nondimensional depth, del0=[lam,-mu], char=[L_ice, T_ice]
%% (L_ice=characteristic length of the thicker ice) & L is the natural length. 

if length(Z)==1
	T=Z; Z=NDphyspram([]); Z{2}=T;
elseif length(Z)==2
	T=Z(1); theta=Z(2);
	Z=NDphyspram([]); Z{2}=T; Z{4}=theta;
elseif length(Z)==3
	T=Z{1}; theta=Z{2}; H_dim=Z{3};
	Z={1,T,H_dim,theta};
end

%% get impt parameters
y=ND1(Z,N);%={gam1,H,th,del0,char,L}
H=y{2}; th=pi*y{3}/180; del0=y{4}; char=y{5}; L=y{6};

%% sort out roots
gam1=y{1};
gam0vec=gam1(1)*[cos(th),sin(th)];
gam1=sign((imag(gam1)>=0)-.5).*gam1;

if nargin==3
	gam1={gam1,alp1,H,el,del0,L,char};
end
