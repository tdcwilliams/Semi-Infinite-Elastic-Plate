function s=GEN_intrinsic_admittance(Z,hr);
%% CALL: s=GEN_intrinsic_admittance(Z,[hr0 hr1]);
%% Z={[lam,-mu],[gam0 gam1],el,[H0 H1]};
%% gam0 sat's (hr0^3*gam0^4+lam-hr0*mu)*gam0*tanh(gam0*H0)=1;
%% gam1 sat's (hr1^3*gam1^4+lam-hr1*mu)*gam1*tanh(gam1*H1)=1.

del0=Z{1}; gamgam=Z{2}; el=Z{3}; HH=Z{4};%gamgam
lam=del0(1); mu=-del0(2);
gam0=gamgam(1); gam1=gamgam(2);
H0=HH(1); H1=HH(2);
%%
hr0=hr(1); hr1=hr(2);
if el>=gam1
	s=0;
else
	alp0=sqrt(gam0^2-el^2);
	alp1=sqrt(gam1^2-el^2);
	Lam0=hr0^3*gam0^4+lam-hr0*mu;
	Lam1=hr1^3*gam1^4+lam-hr1*mu;
%a0=gam0/alp0*GEN_calc_residue({lam,mu,H0},hr0,gam0)
%	a2=gam1/alp1*GEN_calc_residue({lam,mu,H1},hr1,gam1)
	AG0=-gam0/alp0*GEN_calc_residue({lam,mu,H0},hr0,gam0)*Lam0^2;
	AG1=-gam1/alp1*GEN_calc_residue({lam,mu,H1},hr1,gam1)*Lam1^2;
	s=AG0/AG1;
end
