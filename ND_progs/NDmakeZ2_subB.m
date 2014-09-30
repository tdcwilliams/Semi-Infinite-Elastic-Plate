function [Rts,rts,HHsig,el,del0,L]=NDmakeZ2_subB(Z,hh,Nmax)
%% Z=T or [T th] or {T H_dim th}

%% Z->{T,th,H_dim};
if length(Z)==1
	T=Z; Z_=NDphyspram([]);
	Z={T}; Z(2:3)=Z_([4 3]);
elseif length(Z)==2
	T=Z(1); theta=Z(2);
	Z_=NDphyspram([]); Z={T,theta,Z_{3}};
elseif length(Z)==3
	Z(2:3)=Z([3 2]);
end

%% get impt parameters:
[Rts,HHsig,th,del0,L]=ND2subB(Z,hh,Nmax);

%% sort out roots:
gam0=Rts{1}(1); el=gam0*sin(th);
for j=1:3
	gam=Rts{j}; alp=sqrt(gam.^2-el^2);
	rts{j}=sign((imag(alp)>=0)-.5).*alp;
end
