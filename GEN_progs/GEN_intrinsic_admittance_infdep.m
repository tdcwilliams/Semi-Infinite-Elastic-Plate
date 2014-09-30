function s=GEN_intrinsic_admittance_infdep(Z,hr);
%% CALL: s=GEN_intrinsic_admittance(Z,[hr0 hr1]);
%% s=intrinsic admittance for wave going from region 0 to region 1
%% Z={[lam,-mu],[gam0 gam1],el};
%% gam0 sat's (hr0^3*gam0^4+lam-hr0*mu)*gam0=1;
%% gam1 sat's (hr1^3*gam1^4+lam-hr1*mu)*gam1=1.

del0=Z{1};
gamgam=Z{2};
el=Z{3};
lam=del0(1);
mu=-del0(2);
gam0=gamgam(1);
gam1=gamgam(2);
%%
hr0=hr(1);
hr1=hr(2);
if el>=gam1
  s=0;
else
  alp0=sqrt(gam0^2-el^2);
  alp1=sqrt(gam1^2-el^2);
  AG0=-gam0/alp0/(4*hr0^3*gam0^5+1);
  AG1=-gam1/alp1/(4*hr1^3*gam1^5+1);
  s=AG0/AG1;
end
