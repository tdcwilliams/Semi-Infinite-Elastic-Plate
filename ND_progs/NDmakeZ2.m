function [Rts,rts,H,el,del0,L,char]=NDmakeZ2(Z,hr,N)
%% MAIN NONDIMENSIONALISATION PROGRAM WHEN 2 SHEET THICKNESSES ARE INVOLVED.
%% Differs from NDmakeZ2_rel.m in that it nondimensionalises with respect to
%% the lh ice sheet, whereas *rel.m ND's wrt the thicker of the 2 sheets.
%%
%% CALL: [Rts,rts,H,el,del0,L,char]=NDmakeZ2(Z,hr,N)
%%
%% INPUTS: Z=T or [T,th] or {h,T,H_dim} or {h,T,H_dim,th},
%% hr=h1/h0 (hr=0 => open water to far right),
%% N is the number of imaginary roots wanted.
%%
%% OUTPUTS: Rts={K,gam}, where K is the (N+3)-vector of the roots
%% satisfying the dispersion relation for the lhs, gam is the
%% (N+3)- or )N+1)-vector satisfying the dispersion relation for the rhs;
%% el=K(1)*sin(theta);
%% rts={kn,alp}, where kn=sqrt(K.^2-el^2) and alp=sqrt(gam.^2-el^2)
%% are the wavenumbers in the x direction
%% (square roots are taken from the upper complex half-plane);
%% H is the nondimensional depth, del0=[lam,-mu], L is the natural length
%% and char=[L_ice, T_ice], where L_ice is the characteristic length
%% for the lh ice and T_ice is its characteristic time.

if length(Z)==1
	T=Z; Z=NDphyspram([]); Z{2}=T;
elseif length(Z)==2
	T=Z(1); theta=Z(2);
	Z=NDphyspram([]); Z{2}=T; Z{4}=theta;
elseif length(Z)==3
	theta=0; Z{4}=theta;	
end

y=ND2(Z,hr,N);%% y={K,gam,H,th,del0,char,L}
H=y{3}; th=pi*y{4}/180; del0=y{5}; char=y{6}; L=y{7};
gam1=y{1}; gam2=y{2}; Rts={gam1,gam2}; el=gam1(1)*sin(th);
alp1=sqrt(gam1.^2-el^2); alp1=sign((imag(alp1)>=0)-.5).*alp1;
alp2=sqrt(gam2.^2-el^2); alp2=sign((imag(alp2)>=0)-.5).*alp2;
rts={alp1,alp2};
if nargout==1
	Rts={Rts,rts,H,el,del0,L,char};
end
