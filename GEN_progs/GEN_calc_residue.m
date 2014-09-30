function res=GEN_calc_residue(Z,hr,gam);
%% CALL res=GEN_calc_residue(Z,hr,gam);
%% where Z={lam,mu,H}, gam_n is a root of the dispersion relation
%% f(gam)=coth(gam*H)/gam-(Dr*gam^4+lam-hr*mu);
%% y=1/f'(gamma).
lam=Z{1}; mu=Z{2}; H=Z{3}; Dr=hr^3; mr=hr;
Lam=Dr*gam.^4+lam-mr*mu; dLam_times_gam=Lam+4*Dr*gam.^4;
denom=H*(Lam.^2.*gam.^2-1)+dLam_times_gam;
res=-gam./denom;
