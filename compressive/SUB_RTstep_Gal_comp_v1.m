function [R,T,y] = SUB_RTstep_Gal_comp_v1(...
			phys_vars,hh,bc,NN,youngs,Ainc_vec)

%% Notation closer to Williams (2014);
%%
%% calc's scattering coefficients for an ice step
%% CALL: [R,T,y]=...
%%   SUB_RTstep_Gal_comp_v1(phys_vars,hh,bc,NN,INC_SUB)
%% INPUTS: N is no of imag roots to use;
%%   phys_vars=(vector) T or [T, theta_inc=angle of incidence];
%%     or (cell) {T,theta_inc,H_dim}.
%%   where H_dim is distance from sea floor
%%     to bottom of thickest ice sheet;
%%   NN=vector[Nterms,Nroots]
%%     -> no of polys to use for P' in Galerkin expansion
%%     & no of roots to use
%% OUTPUTS: R&T are reflection and transmission coefficients;
%%   y=...
%%TODO Young's -> EE

INC_SUB  = 1;
do_test  = 1;
if ~exist('phys_vars')
   period      = 10;%% wave period [s]
   theta_inc   = 0;%% wave incident angle [degrees]
   H_dim       = 100;%% water depth [m]
   phys_vars   = {period,theta_inc,H_dim};
end
if ~exist('hh')
   hh = [1 2];
   %hh = [2 1];
   %hh = [0 1];
   %hh = [1 0];
   %hh = [1 1];
end
if ~exist('bc')
   %bc = 1;%% free edge conditions;
   bc = 0;%% frozen edge conditions;
end
if ~exist('NN')
   NN = [10 1000];%% [N_poly, N_roots];
end
if ~exist('youngs')
   youngs   = 5.45e9*[1 1];
   %youngs   = 5.45e9*[1 .8];
   %youngs   = [1 1]*NDphyspram(1);
elseif length(youngs)==1
   youngs   = [1 1]*youngs;
end

%%incident waves:
if ~exist('Ainc_vec')
   Ainc2    = 1;%%f-g wave from lhs;
   Binc2    = 0;%%f-g wave from rhs;
   Ac_inc   = 0;%%compressive wave from lhs;
   Bc_inc   = 0;%%compressive wave from lhs;
   Ainc_vec = [Ac_inc,Ainc2,Bc_inc,Binc2];
end
Ac_inc   = Ainc_vec(1);
Ainc2    = Ainc_vec(2);
Bc_inc   = Ainc_vec(3);
Binc2    = Ainc_vec(4);

if nargin==0
   %% do some tests
   %% - print |R|&|T|, and check energy;
   do_test  = 1;
end

%% tried to improve the number of roots needed
%% but not working at the moment - so don't recommend this option;
FAST_KERNEL = 0;
if exist('Mlog')
   FAST_KERNEL = 1;
end

%% DO NONDIMENSIONALIZATION:
Nterms   = NN(1);
if length(NN)==1
  decay_exp = 2;
  ndp       = 8;
  Href      = [];%%automatic Href is H_dim
  Ninput    = {decay_exp,ndp,Href};
    %% inputting info about Nroots in this way lets the
    %% nondim program determine how many are needed for a
    %% series of O[(n*pi*Href/H)^(-decay_exp)]
    %% to converge to ndp decimal places
else
  Nroots = NN(2);
  Ninput = Nroots*ones(1,1+INC_SUB);
end

%[gam,alp,hw,alpy,L]=NDmakeZ_wtr(phys_vars,Nroots);
if INC_SUB==1
  [Rts,rts,HH,alpy,del0,L] = NDmakeZ2_sub(phys_vars,hh,Ninput,youngs);
else
  [Rts,rts,H,alpy,del0,L]  = NDmakeZ2_rel(phys_vars,hh,Ninput);
  HH                       = [H H];
end

%%want larger ice thickness on right:
DO_SWAP  = ( hh(1)>hh(2) );
if DO_SWAP
   hh       = fliplr(hh);
   youngs   = fliplr(youngs);
   Rts      = Rts([2 1]);
   rts      = rts([2 1]);
   HH       = HH([2 1]);
   %%
   tmp   = Ainc2;
   Ainc2 = Binc2;
   Binc2 = tmp;
   %%
   tmp      = Ac_inc;
   Ac_inc   = Bc_inc;
   Bc_inc   = tmp;
end

if ~DO_SWAP
   Z_out = {{Rts{1}/L,Rts{2}/L},{rts{1}/L,rts{2}/L},...
             HH*L,alpy/L,del0,L};
else
   Z_out = {{Rts{2}/L,Rts{1}/L},{rts{2}/L,rts{1}/L},...
             HH([2 1])*L,alpy/L,del0,L};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Elastic constants:
nu       = NDphyspram(5);
E1       = youngs(1);%%E the same on both sides
E2       = youngs(2);
rho_wtr  = NDphyspram(3);
rho_ice  = NDphyspram(4);
om       = 2*pi/phys_vars{1};

%%Compressional stuff for u problem:
mu1_lame = E1/2/(1+nu);
mu2_lame = E2/2/(1+nu);
Kc1_dim  = E1*hh(1)/(1-nu^2);%%compressional rigidity ~ Pa*m ~ rho_wtr*om^2*L^2*L
Kc2_dim  = E2*hh(2)/(1-nu^2);
Kc1      = Kc1_dim/(rho_wtr*om^2*L^2)/L;
Kc2      = Kc2_dim/(rho_wtr*om^2*L^2)/L;
%%
m1_dim   = rho_ice*hh(1);
m2_dim   = rho_ice*hh(2);
sig1_dim = m1_dim/rho_wtr;
sig2_dim = m2_dim/rho_wtr;
zc1_dim  = -sig1_dim+hh(1)/2;
zc2_dim  = -sig2_dim+hh(2)/2;
dzc_dim  = zc1_dim-zc2_dim;%pause
%%
sig1  = sig1_dim/L;
sig2  = sig2_dim/L;
zc1   = zc1_dim/L;
zc2   = zc2_dim/L;
dzc   = dzc_dim/L;
%%
if Kc1==0
   kc1_dim  = 0;
else
   kc1_dim  = sqrt(om^2*m1_dim/Kc1_dim);
end
kc2_dim  = sqrt(om^2*m2_dim/Kc2_dim);
kc1      = kc1_dim*L;
kc2      = kc2_dim*L;
%%
M1_rel   = rho_wtr*om^2*L^4;%%scale bending moment by this, since \int\sig_ij.(z-z_c).dz ~ Pa*m^2 ~ rho_wtr*om^2*L^2*L^2 = D1/L
M0_rel   = rho_wtr*om^2*L^3;%%scale zero-th moment by this, since \int\sig_ij.dz ~ Pa*m ~ rho_wtr*om^2*L^2*L = D1/L^2
J1_dim   = rho_ice*hh(1)^3/12;
J2_dim   = rho_ice*hh(2)^3/12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
hr = hh/hh(2);
Dr = [E1/E2,1].*hr.^3;
%%
gam1  = Rts{1};
alp1  = rts{1};
H1    = HH(1);
gam2  = Rts{2};
alp2  = rts{2};
H2    = HH(2);

%%make real parts of gam be >0 (for M0_wtr,M1_wtr)
jneg        = find(real(gam1)<0);
gam1(jneg)  = -gam1(jneg);
jneg        = find(real(gam2)<0);
gam2(jneg)  = -gam2(jneg);

if 0%%test rts etc
  HH
  for j=1:2
    gg            = Rts{j};
    tst_disprel   = ( Dr(j)*gg.^4+del0*[1;hr(j)] ).*...
                     gg.*tanh(gg*HH(j))-1
  end
  return
end

%%
lam   = del0(1);
mu    = -del0(2);
nu1   = (1-nu)*alpy^2;
sig2  = mu;
H     = H2+sig2;
%%
gam0  = gam1(1);
alp0  = alp1(1);
%%
BGzz1 = calc_res({lam,mu,H1},Dr(1),hr(1),gam1).*gam1./alp1;%%BG_0^(0)
Lam1  = Dr(1)*gam1.^4+lam-hr(1)*mu;
BGz1  = -Lam1.*BGzz1;%%BG_0^(1); why '-'? (in W&P paper, but seems wrong)
BGz1B = Lam1.*BGzz1;%%fix minus sign in BG_0^(1)
BG1   = Lam1.*BGz1B;%%BG_0^(2)
%%
BGzz2 = calc_res({lam,mu,H2},Dr(2),hr(2),gam2).*gam2./alp2;
Lam2  = Dr(2)*gam2.^4+lam-hr(2)*mu;
BGz2  = -Lam2.*BGzz2;%%why '-'? (in W&P paper, but seems wrong)
BGz2B = Lam2.*BGzz2;%%fix minus sign
BG2   = Lam2.*BGz2B;
%%

%% CALC Fm0 & Fmr:
mvec  = (0:Nterms)';
alpC  = .5-1/3*INC_SUB*(hr(1)~=1);%%=1/6 if submergence included;
			       %%=1/2 if no submergence.
kap1     = -1i*gam1*H2;
besJ1    = besselj(2*mvec'+alpC,kap1).';
c_left   = gamma(alpC)*(alpC+2*mvec).*(-1).^mvec;
c_rt1    = (2./kap1).^alpC./cosh(gam1*H1);
F1       = diag(c_left)*besJ1*diag(c_rt1);%%same as (37)
%%
kap2  = -1i*gam2*H2;
besJ2 = besselj(2*mvec'+alpC,kap2).';
c_rt2 = (2./kap2).^alpC./cos(kap2);
F2    = diag(c_left)*besJ2*diag(c_rt2);
%%
if 0%%test F1:
   j            = 13;
   np           = 100;
   Np           = 1+2*Nterms;
   tt           = (0:np)'/np;
   varf         = cosh(gam1(j)*H2*tt)/cosh(gam1(j)*H1);
   an           = zeros(1,Np);
   an(1:2:Np)   = F1(:,j);
   varf_ap      = OP_interp_gegenbauer(tt,alpC,an);
   plot(tt,[real(varf),imag(varf)]), hold on;%pause
   plot(tt,[real(varf_ap),imag(varf_ap)],'--r'), hold off;
   return
elseif 0%%test F2:
   j            = 20;
   np           = 100;
   Np           = 1+2*Nterms;
   tt           = (0:np)'/np;
   varf         = cosh(gam2(j)*H2*tt)/cosh(gam2(j)*H2);
   an           = zeros(1,Np);
   an(1:2:Np)   = F2(:,j);
   varf_ap      = OP_interp_gegenbauer(tt,alpC,an);
   plot(tt,[real(varf),imag(varf)]), hold on;%pause
   plot(tt,[real(varf_ap),imag(varf_ap)],'--r'), hold off;
   return
end

if 0
   np           = 100;
   Np           = 1+2*Nterms;
   tt           = (0:np)'/np;
   uu           = alp1(1)*cosh(gam1(1)*H2*tt)/cosh(gam1(1)*H1);
   um           = alp1(1)*F1(:,1);
   an           = zeros(1,Np);
   an(1:2:Np)   = um;
   plot(tt,uu), hold on;
   plot(tt,GEN_split_ri(OP_interp_gegenbauer(tt,alpC,an)),'--r');
   hold off;
   return;
end

if 0%%check iprules
   jt=1:6;
   if 1
     prams1  = {[Dr(1) Dr(1)],lam,-del0(2)*hr(1)*[1 1],H1};
     M1      = GEN_iprules_ice(gam1,gam1,prams1);
     tst1    = M1(jt,jt)
     M1b     = GEN_iprules_ice_num(gam1,gam1,[0 H1],[H1,H1],[H1 H1]);
     tst1b   = M1b(jt,jt)
   else
     prams2  = {[Dr(2) Dr(2)],lam,-del0(2)*hr(2)*[1 1],H2};
     M2      = GEN_iprules_ice(gam2,gam2,prams2);
     tst2    = M2(jt,jt)
     M2b     = GEN_iprules_ice_num(gam2,gam2,[0 H2],[H2,H2],[H2 H2]);
     tst2b   = M2b(jt,jt)
   end
   return;
end
  
%%MAIN KERNEL MATRIX:
if 1%~FAST_KERNEL%%basic way:-other ways don't work
   MK = F2*diag(1i*BG2)*F2.'+F1*diag(1i*BG1)*F1.';%%looks like same as (38c)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Moments of pressure M0w and M1w:
%%\int_{-\sig2}^{-sig1}\varf.dz
denom    = 1+exp(-2*gam1*H1);
num0     = 1-exp(-2*gam1*H2);
exfac    = exp(-gam1*(H1-H2));
M0_wtr   = 1./gam1.^2./Lam1+...
            -1./gam1.*exfac.*num0./denom;

if 0
   tstM0    = tanh(gam1*H1)./gam1-1./gam1.*exfac.*num0./denom;
   tst_dr1  = [gam1.*tanh(gam1*H1),1./Lam1];{tst_dr1}
   tst_dr2  = [gam2.*tanh(gam2*H2),1./Lam2];
   tst_dr1(1:5,:)
   tst_dr2(1:5,:)
   tstM0(1:5)
   M0_wtr(1:5),dzc,H1,H2
   disp('check new dispersion relation progs! (E1,E2)')
   return;
end
%%
%%\int_{-\sig2}^{-sig1}(z-z_c)\varf.dz, z=z_c is center of rhs plate
M1a      = (H1./Lam1-1)./gam1.^2;
num1     = (1-gam1*H2)+(1+gam1*H2).*exp(-2*gam1*H2);
M1b      = num1.*exfac./gam1.^2./denom;
M1_wtr   = M1a+M1b-(H+zc2)*M0_wtr;

if 0%%test M0_wtr,M1_wtr
   jt    = 1:10;
   x1    = gam1*H1;
   x2    = gam1*H2;
   int0  = (sinh(x1)-sinh(x2))./gam1./cosh(gam1*H1);
   int1a = ( x1.*sinh(x1)-cosh(x1) )./gam1.^2./cosh(x1);
   int1b = ( -x2.*sinh(x2)+cosh(x2) )./gam1.^2./cosh(x1);

   %int1_0   = int1-H1*int0;
   %[M0_wtr(jt),int0(jt)]
   %[M1_wtr_0(jt),int1_0(jt)]

   [int1a(jt),M1a(jt)]
   [int1b(jt),M1b(jt)]
   %%
   jt2   = 3;
   f0_in = inline('cosh(k.*(z+H))./cosh(k.*(H-sig1))');%%f0(H,k,sig1,z)
   q0    = quad(@(z)f0_in(H,alp1(jt2),sig1,z),-sig2,-sig1);
   f1_in = inline('(z-zc2).*cosh(k.*(z+H))./cosh(k.*(H-sig))');%%f1(H,k,sig1,z,zc2)
   q1    = quad(@(z)f1_in(H,alp1(jt2),sig1,z,zc2),-sig2,-sig1);
   [q0,M0_wtr(jt2)]
   [q1,M1_wtr(jt2)]

   %%test quad :)
   %f2_in = inline('cosh(a.*b.*c.*d.*z)');
   %q2    = quad(@(z)f2_in(1,1,1,1,z),0,1);
   %f3_in = inline('2*sinh(e.*z).*cosh(a.*b.*c.*d.*z)');
   %q3    = quad(@(z)f3_in(1,1,1,1,1,z),0,1);
   %[q2,sinh(1)]
   %[q3,(cosh(2)-1)/2]
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%contribution from side integrals to rr2;
%%Q2 (psi(0+)) column = #5:
rr_wtr_psi2 = -2i*BG1.*M1_wtr;
%%U2 (u(0+)) column = #7:
rr_wtr_U2   = -2i*BG1.*M0_wtr;
jside       = [5 7];

%%matrix M0w=M_M0w*rr2
M_M0w = -M0_wtr.';%%\int\sig_11 dz=\int(-P)dz=-\int\phi dz
M_M1w = -M1_wtr.';%%\int\sig_11*(z-zc1)dz=\int(-P)*(z-zc1)dz=-\int\phi*(z-zc1)dz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%FORCING TERMS:
fm1   =  gam1.^2-nu1;
fm2   =  gam2.^2-nu1;
E1    =  [1+0*alp1,-Dr(1)*fm1];
E2    =  [1+0*alp2,-Dr(2)*fm2];
E1B   =  [-1+0*alp1,Dr(1)*fm1];%% S=-(L^2/D1)*D*w_xxx=D/D1*\pa_{xND}^3(wND), psi=-w_x
E2B   =  [-1+0*alp2,Dr(2)*fm2];



%%SOLVE INTEGRAL EQN:
%Ainc  = -1i;
Ainc  = 1;
ME2   = F2*diag(1i*BGz2)*E2;
ME1   = F1*diag(1i*BGz1)*E1;
uu    = MK\[Ainc*F1(:,1),ME1,ME2];%%(38a)
ME2B  = F2*diag(1i*BGz2B)*E2B;
ME1B  = F1*diag(1i*BGz1B)*E1B;

%%extra forcing from side integral terms in rr2;
forcing_u            = [Ainc2*F1(:,1)-Binc2*F2(:,1),ME1B,ME2B,0*ME1B];
forcing_u(:,jside)   = forcing_u(:,jside)+F1*[rr_wtr_psi2,rr_wtr_U2];
uuB                  = MK\forcing_u;
%%
rr          = 2*diag(BG1)*F1.'*uu;%%in (24) of W&P: rr=1i*R_0(\alpha)
rr(1)       = rr(1)+1i*Ainc;
rr(:,2:3)   = rr(:,2:3)-2*diag(BGz1)*E1;
%%
rr2         = -2i*diag(BG1)*F1.'*uuB;%%now same as (24) of W&P (but Q vector -> -Q)
rr2(1)      = rr2(1)+Ainc2;%%now same as (24) of W&P
rr2(:,2:3)  = rr2(:,2:3)+2i*diag(BGz1B)*E1B;

%%from side integrals;
rr2(:,jside)   = rr2(:,jside)+[rr_wtr_psi2,rr_wtr_U2];

%%
tt          = -2*diag(BG2)*F2.'*uu;
tt(:,4:5)   = tt(:,4:5)+2*diag(BGz2)*E2;
tt2         = 2i*diag(BG2)*F2.'*uuB;
tt2(1)      = tt2(1)+Binc2;
tt2(:,4:5)  = tt2(:,4:5)-2i*diag(BGz2B)*E2B;%%like W&P paper (but Q vector -> -Q)
if 0
   jt = 1:5;
   rr(jt,:)/1i, rr2(jt,:)
   tt(jt,:)/1i, tt2(jt,:),return
end
%%
M_PM1 = E1.'*diag(-1./Lam1);%%why '-'? (in W&P paper, but seems wrong)
M_PM2 = E2.'*diag(-1./Lam2);%%why '-'? (in W&P paper, but seems wrong)
M_PM1B   = diag([-1,1])*E1B.'*diag(1./Lam1);%% M=-Dr*Lm(w)=-Dr*[-fm*(wND)]=Dr*fm
M_PM2B   = diag([-1,1])*E2B.'*diag(1./Lam2);%%

%%APPLY EDGE CONDITIONS:
if Dr(1)==0%%water on left:
   j_dis    = 2:4
   j_dis2   = [2 3 4 6];%%no more [S1 (2),Q1 (3),U1 (6)], & S2=0 (4)
                            %%keep Q2 (5) and U2 (7)
   rr(:,j_dis)    = [];
   tt(:,j_dis)    = [];
   rr2(:,j_dis2)  = [];
   tt2(:,j_dis2)  = [];
   
   %%old
   M2  = M_PM2(2,:)*tt;
   Q2  = -M2(1)/M2(2);
   
   %%bending moment eqn
   M2B      = M_PM2B(2,:)*tt2;
   M2B(1)   = M2B(1)+Binc2*M_PM2B(2,1);%%'+' since M has even derivatives
   M2Bw     = M_M1w*rr2;
   M2Bw(1)  = M2Bw(1)+M_M1w(1)*Ainc2;
   Medge    = M2B-M2Bw

   %%compression at edge: U2=(Bc_inc+Bc_scat)
   M0w      = M_M0w*rr2;
   M0w(1)   = M0w(1)+ M_M0w(1)*Ainc2;
   %%
   Medge(2,:)  = M0w/(1i*kc2*Kc2);
   Medge(2,1)  = Medge(2,1)+2*Bc_inc;
   Medge(2,3)  = Medge(2,3)-1;



   Q2B      = -Medge(:,2:3)\Medge(:,1);%pause
   %%
   rn  = rr*[1;Q2];
   tn  = tt*[1;Q2];
   rn2      = rr2*[1;Q2B];
   tn2      = tt2*[1;Q2B];
   Bc_scat  = Q2B(2)-Bc_inc;
   %%
   Ac_inc   = 0;
   Ac_scat  = 0;
   %%
   v_unk       = [1;Q2];
   v_unk2      = eye(4,1);
   v_unk2(3:4) = Q2B;%%Q2,U2
   uu(:,j_dis) = [];
elseif bc==1%%free edges: col's of Medge corresp to [1 psi1 psi2 U1 U2]
   rr(:,[2 4]) = [];
   tt(:,[2 4]) = [];%% S(L)=S(R)=0
   rr2(:,[2 4])   = [];
   tt2(:,[2 4])   = [];%% S(L)=S(R)=0
   %%
   MM    = [M_PM1(2,:)*rr; M_PM2(2,:)*tt];
   MM(1) = MM(1)+1i*Ainc*M_PM1(2,1);
   QQ    = -MM(:,2:3)\MM(:,1);%%[w_x L,w_x R]
   
   %%Bending moment eqn (LHS)
   Medge       = M_PM1B(2,:)*rr2;%%bending moment: [L;R]
   Medge(1,1)  = Medge(1,1)+Ainc2*M_PM1B(2,1);%%need to add the LHS incident wave bending moment

   %%bending moment eqn (RHS)
   M2B         = M_PM2B(2,:)*tt2;
   M2B(1)      = M2B(1)+Binc2*M_PM2B(2,1);%%'+' since M has even derivatives
   M2Bw        = M_M1w*rr2;
   M2Bw(1)     = M2Bw(1)+M_M1w(1)*Ainc2;
   Medge(2,:)  = M2B-M2Bw;

   %%compression at edge (LHS):
   %% *U1=Ac_inc+Ac_scat=2*Ac_inc, or U1-2*Ac_inc=0
   Ac_scat     = Ac_inc;%%u_x=0
   U1          = 2*Ac_inc;
   Medge(3,1)  = U1;
   Medge(3,4)  = -1;

   %%compression at edge (RHS): Bc_scat-Bc_inc = U2-2*Bc_inc = M0w/(1i*kc2*Kc2)
   M0w      = M_M0w*rr2;
   M0w(1)   = M0w(1)+ M_M0w(1)*Ainc2
   M0i      = 1i*kc2*Kc2*[-2*Bc_inc,0,0,0,1];
   %%
   Medge(4,:)  = M0w-M0i
   %pause
   %Medge(4,:)/(1i*kc2*Kc2),return
   %%
   v_unk          = eye(5,1);
   v_unk([3 5])   = QQ;
   rn             = rr*[1;QQ];
   tn             = tt*[1;QQ];
   %%
   QQ2               = -Medge(:,2:end)\Medge(:,1);%%[psi1 psi2 U1 U2]
   v_unk2            = eye(7,1);
   v_unk2([3 5 6 7]) = QQ2;
   rn2               = rr2*[1;QQ2];
   tn2               = tt2*[1;QQ2];
   U2                = QQ2(4);
   Bc_scat           = U2-Bc_inc;
   M0w*[1;QQ2]
   Bc = {Bc_inc,Bc_scat},pause

%  {M0w*[1;QQ2],...
%     M0i*[1;QQ2],...
%     M_M0w*rn2+ M_M0w(1)*Ainc2,...
%     1i*kc2*Kc2*(Bc_scat-Bc_inc)},pause

   if 0
      %[M_PM1(2,:)*rr; M_PM2(2,:)*tt]/1i
      %[M_PM1(2,:)*rr2; M_PM2(2,:)*tt2]
      MM/1i,MM2
      QQ,QQ2,return
   end
else%%frozen edges:[1,S1,Q1,S2,Q2,U1,U2]->[1,S1,Q1,U1,U2]
   rr(:,2:3)   = rr(:,2:3)+rr(:,4:5);%% [-S,w_x](L)=[-S,w_x](R) so eliminate (R)
   rr(:,4:5)   = [];


   %%
   tt(:,2:3)   = tt(:,2:3)+tt(:,4:5);
   tt(:,4:5)   = [];
   %%
   PM1      = M_PM1*rr;
   PM1(:,1) = PM1(:,1)+1i*Ainc*M_PM1(:,1);
   PM2      = M_PM2*tt;
   dPM      = PM2-PM1;
   Q1       = -dPM(:,2:3)\dPM(:,1);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%S2=S1, psi2=psi1 (cty of horizontal displacement pt1)
   rr2(:,2:3)  = rr2(:,2:3)+rr2(:,4:5);%% [-S,w_x](L)=[-S,w_x](R) so eliminate (R)
   rr2(:,4:5)  = [];
   tt2(:,2:3)  = tt2(:,2:3)+tt2(:,4:5);
   tt2(:,4:5)  = [];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% cty of horizontal displacement pt2
   Medge = [0 0 dzc -1 1];
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% 0 moment (RHS): (1i*kc2*Kc2)*(Bc_scat-Bc_inc) = (1i*kc2*Kc2)*(U2-2*Bc_inc)
   M0i2     = 1i*kc2*Kc2*[-2*Bc_inc,0,0,0,1];

   %% 0 moment (LHS): (1i*kc1*Kc1)*(Ac_inc-Ac_scat) = (1i*kc1*Kc1)*(2*Ac_inc-U1)
   M0i1     = 1i*kc1*Kc1*[2*Ac_inc,0,0,-1,0];

   %% 0 moment from water:
   M0w      = M_M0w*rr2;
   M0w(1)   = M0w(1)+ M_M0w(1)*Ainc2;

   %% full 0 moment eqn:
   Medge(2,:)  = M0w+M0i1-M0i2;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Bending moment (LHS)
   M1i1     = M_PM1B(2,:)*rr2;%%bending moment: [L;R]
   M1i1(1)  = M1i1(1)+Ainc2*M_PM1B(2,1);%%need to add the LHS incident wave bending moment

   %% Bending moment (RHS)
   M1i2     = M_PM2B(2,:)*tt2;
   M1i2(1)  = M1i2(1)+Binc2*M_PM2B(2,1);%%'+' since M has even derivatives

   %% Bending moment (water)
   M2Bw        = M_M1w*rr2;
   M2Bw(1)     = M2Bw(1)+M_M1w(1)*Ainc2;
   Medge(3,:)  = M1i2-M1i1-M2Bw-dzc*M0i1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% w (LHS):
   M_W1        = M_PM1B(1,:)*rr2;
   M_W1(1)     = M_W1(1)+Ainc2*M_PM1B(1,1);

   %% w (RHS):
   M_W2        = M_PM2B(1,:)*tt2;
   M_W2(1)     = M_W2(1)+Binc2*M_PM2B(1,1);

   %% w cts:
   Medge(4,:) = M_W2-M_W1
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   QQ2   = -Medge(:,2:end)\Medge(:,1);
   S1    = QQ2(1);
   psi1  = QQ2(2);
   U1    = QQ2(3);
   U2    = QQ2(4);
   %%
   v_unk2   = [1;S1;psi1;S1;psi1;U1;U2];
   rn2      = rr2*[1;QQ2];
   tn2      = tt2*[1;QQ2];
   Ac_scat  = U1-Ac_inc;
   Bc_scat  = U2-Bc_inc;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %% old
   v_unk = [1;Q1;Q1];
   rn    = rr*[1;Q1];
   tn    = tt*[1;Q1];

   tstM0 = -Medge(2,:)*[1;QQ2]
   {M0w*[1;QQ2],M0i1*[1;QQ2],M0i2*[1;QQ2]}

   if 0
      %[M_PM1(2,:)*rr; M_PM2(2,:)*tt]/1i
      %[M_PM1(2,:)*rr2; M_PM2(2,:)*tt2]
      %MM/1i,MM2
      Q1,Q1B,return
   end
end

if 0
   jt=1:10;
   [v_unk,v_unk2]
   [rn(jt)/1i,rn2(jt)]
   [tn(jt)/1i,tn2(jt)],return
end

%%REFLECTION AND TRANSMISSION COEFFICIENTS:
if DO_SWAP==0
   R2 = rn2(1)/Ainc2;
   T2 = tn2(1)/Ainc2;
   s  = ( BG1(1)/BG2(1) );%intrinsic admittance
   y2 = {Z_out,rn,tn,s};
else
   R2 = tn2(1)/Binc2;
   T2 = rn2(1)/Binc2;
   s  = ( BG2(1)/BG1(1) );%intrinsic admittance
   y2 = {Z_out,tn,rn,s};
   %%
   %v_unk2   = v_unk2([1 4 5 2 3]);
end
R  = -1i*rn(1)/Ainc;
T  = -1i*tn(1)/Ainc;
s  = ( BG1(1)/BG2(1) );%intrinsic admittance
y  = {Z_out,rn,tn,s};

%if DO_SWAP==0%%test R,T for LHS wave
%   [R,R2]
%   [T,T2]
%   return
%else
%   Rm  = -R'*T/T';
%   Tm  = (1-abs(R)^2)/T';
%   [Rm,R2]
%   [Tm,T2]
%   return
%end

if DO_SWAP
   R   = -R'*T/T';
   T   = (1-abs(R)^2)/T';
   y(2:3) = [];%%can't just swap, but there is a way - FIX!
   %%
   %y(2:3) = y([3 2]);%%swap rn,tn
end

if 0%%full Green's fxn formulation:
   %%vm=\phi(0-,z)=\sum r_n\varf_n(z)
   uu2   = uu*v_unk;%%u(z) coefficients
   S0    = -v_unk(1+1);%%in W&P paper, S is defined as -1*[shear resultant]
   psi0  = -v_unk(1+2);%%\psi=-w_x on left 

   W_Gz1    = 1i*BGz1;%%G_z
   psi_Gz1  = 1i*(-1i*alp1).*BGz1;%%-G_zx
   M_Gz1    = 1i*Dr(1)*alp1.^2.*BGz1;%%-D0*G_zxx
   S_Gz1    = -1i*Dr(1)*(1i*alp1).^3.*BGz1;%%-D0*G_zxxx
   %%
   Ainc  = -1i;
   vn0   = Ainc*eye(length(alp1),1)+...
           -2i*diag(BG1)*F1.'*uu2+...
           +2*psi0*M_Gz1+...
           -2*S0*W_Gz1;%%this is rn/1i
   jt = 1:10;
   %[vn0(jt),-1i*rn(jt)]

   PM_1  = M_PM1*rn;
   W0    = PM_1(1); 
   M0    = -PM_1(2); 
   vn1   = vn0+2*W0*S_Gz1+...
            -2*M0*psi_Gz1;%%these are the expansion coefficients for phi(0,z) from Green's thm


   Imn1  = GEN_iprules_ice_num(gam1,gam1,[0 H1],[H1,H1],[H1 H1]);%%explicit integration using cosh product rule
   rn1   = -1i*diag(BG1)*F1.'*uu2+...
            +1i*diag(1i*alp1.*BG1)*Imn1*vn1+...
            +psi0*M_Gz1+...
            -S0*W_Gz1+...
            +2*W0*S_Gz1+...
            -2*M0*psi_Gz1;
   [-1i*rn(jt),vn0(jt),rn1(jt)]


   Imn2  = GEN_iprules_ice_num(gam2,gam2,[0 H2],[H2,H2],[H2 H2]);

   return
end



if do_test==1%%test new results
   disp('R&T, |R|&|T|, |R|^2+s*|T|^2:');
   disp([R2,T2])
   disp([abs(R2),abs(T2)])
   disp(R2*R2'+s*T2*T2')
   %%
   fp1 = gam1.^2+nu1;
   fp2 = gam2.^2+nu1;

   %% New energy check:
   Ew_fac1  = -1/2/BG1(1);
   Ew_fac2  = -1/2/BG2(1);
   %Ew_fac1  = -alp1(1)/2/BG1(1);
   %Ew_fac2  = -alp2(1)/2/BG2(1);
   Ec_fac1  = kc1*Kc1;
   Ec_fac2  = kc2*Kc2;
   %%
   Ec_in    = Ec_fac1*abs(Ac_inc )^2 + Ec_fac2*abs(Bc_inc )^2
   Ec_out   = Ec_fac1*abs(Ac_scat)^2 + Ec_fac2*abs(Bc_scat)^2
   Ew_in    = Ew_fac1*abs(Ainc2  )^2 + Ew_fac2*abs(Binc2  )^2
   Ew_out   = Ew_fac1*abs(rn2(1) )^2 + Ew_fac2*abs(tn2(1) )^2
   tstE_new = [Ew_in+Ec_in,Ew_out+Ec_out]
   %%
   [Ainc2,tn2(1)],[Binc2,rn2(1)]
   if bc==0
      Abc   = {Ac_inc,Bc_scat}
      Bac   = {Bc_inc,Ac_scat}
   else
      dzc
      Ac = {Ac_inc,Ac_scat}
      Bc = {Bc_inc,Bc_scat}
   end
   %return
   GEN_pause;



   %%test definitions of unknown constants;
   an    = rn2./Lam1;
   bn    = tn2./Lam2;
   an0   = Ainc2/Lam1(1);
   bn0   = Binc2/Lam2(1);
   %%
   w_rhs    = bn0+sum(bn);
   w_lhs    = an0+sum(an);
   M_rhs    = -Dr(2)*( (-fm2(1)*bn0)+sum(-fm2.*bn) );%%M=-D*Lm*w -> -D*-fm
   M_lhs    = -Dr(1)*( (-fm1(1)*an0)+sum(-fm1.*an) );
   %return
   %%
   psi_rhs  = -( (-1i*alp2(1)*bn0)+sum(1i*alp2.*bn) );%%-w_x
   psi_lhs  = -( ( 1i*alp1(1)*an0)-sum(1i*alp1.*an) );
   S_rhs    = -Dr(2)*( (-1i*alp2(1).*(-fp2(1))*bn0)+sum(1i*alp2.*(-fp2).*bn) );%%S=-D*Lp*w_x -> -D*(1i*alp)*(-fp)=D*(1i*alp)*(fp)
   S_lhs    = -Dr(1)*( ( 1i*alp1(1).*(-fp1(1))*an0)-sum(1i*alp1.*(-fp1).*an) );
   %return
   wpms1 = [w_lhs;psi_lhs;M_lhs;S_lhs];
   wpms2 = [w_rhs;psi_rhs;M_rhs;S_rhs];

   Mw = [M_M0w*rn2+M_M0w(1)*Ainc2;...
         M_M1w*rn2+M_M1w(1)*Ainc2]

   M0_i1 = 1i*kc1*Kc1*(Ac_inc-Ac_scat)
   M0_i2 = 1i*kc2*Kc2*(Bc_scat-Bc_inc)
   U_1   = Ac_inc+Ac_scat;
   U_2   = Bc_inc+Bc_scat;

   Bc = {Bc_inc,Bc_scat}
   Vunk_chk = [S_lhs;psi_lhs;S_rhs;psi_rhs;U_1;U_2];

   nsp      = 1;
   nsub     = 1;
   DO_PLOT  = 1;
   if DO_PLOT
      if bc==0 & hh(1)>0
         nsp   = nsp+1;
         subplot(2,2,2);
         zp   = linspace(-sig1,hh(1)/L-sig1,40)';
         y1    = real(U_1+psi_lhs*(zp-zc1));
         y2    = real(U_2+psi_rhs*(zp-zc2));
         plot(y1,zp,'k');
         hold on;
         plot(y2,zp,'--r');
         hold off;
         GEN_proc_fig('u/L','z/L');
      end
      
      %% submerged part of ice: phi_x=U_2+(z-zc2)*psi_rhs
      subplot(2,2,1);
      zp2   = linspace(-sig2,-sig1,40)';
      vf2   = cosh((zp2+H)*gam1.')*diag(1./cosh(gam1*H1));
      y21   = real(1i*alp1(1)*Ainc2*vf2(:,1)+...
                  +vf2*(-1i*alp1.*rn2));
      y22   = real(U_2+psi_rhs*(zp2-zc2));
      plot(y21,zp2,'k');
      hold on;
      plot(y22,zp2,'--r');
      hold off;
      GEN_proc_fig('u/L, \phi_x/L','z/L');

      %% phi matching in water:
      subplot(2,2,3);
      zp3   = linspace(-H,-sig2,140)';
      vf3L  = cosh((zp3+H)*gam1.')*diag(1./cosh(gam1*H1));
      vf3R  = cosh((zp3+H)*gam2.')*diag(1./cosh(gam2*H2));
      y31   = real(Ainc2*vf3L(:,1)+...
                  +vf3L*(rn2));
      y32   = real(Binc2*vf3R(:,1)+...
                  +vf3R*(tn2));
      plot(y31,zp3,'k');
      hold on;
      plot(y32,zp3,'--r');
      hold off;
      GEN_proc_fig('\phi/L^2','z/L');

      %% phi_x matching in water:
      subplot(2,2,4);
      y41   = real((1i*alp1(1)*Ainc2)*vf3L(:,1)+...
                  +vf3L*(-1i*alp1.*rn2));
      y42   = real((-1i*alp2(1)*Binc2)*vf3R(:,1)+...
                  +vf3R*(1i*alp2.*tn2));
      plot(y41,zp3,'k');
      hold on;
      plot(y42,zp3,'--r');
      hold off;
      GEN_proc_fig('\phi_x/L','z/L');
   end

   if hh(1)==0
      disp('test definitions of unknown constants')
      [v_unk2(2:4),Vunk_chk([3 4 6])]
      disp('check free edge conditions');
      tstedge    = [ [wpms2(3:4);M0_i2],...
                     [Mw(2);0;Mw(1)]]
   elseif bc==1
      disp('test definitions of unknown constants')
      [v_unk2(2:7),Vunk_chk]
      disp('check free edge conditions');
      tstedge    = [ [wpms1(3:4);wpms2(3:4);M0_i1;M0_i2],...
                     [0;0;Mw(2);0;0;Mw(1)]]
   else
      disp('test definitions of unknown constants')
      [v_unk2(2:7),Vunk_chk]
      disp('check frozen edge conditions');
      tstedge    = [[wpms2-wpms1;M0_i2-M0_i1;U_2-U_1+dzc*psi_rhs],...
                    [0;0;Mw(2)+dzc*M0_i1;0;Mw(1);0]]
      
   end
elseif do_test==2%%test old results

   disp('R&T, |R|&|T|, |R|^2+s*|T|^2:');
   disp([R,T])
   disp([abs(R),abs(T)])
   disp(R*R'+s*T*T')
   %%
   fp1 = gam1.^2+nu1;
   fp2 = gam2.^2+nu1;
   %%
   if bc==1 | hh(1)==0
      disp('check free edge conditions');
      an         = -rn./Lam1;
      bn         = -tn./Lam2;
      Inc        = -1/Lam1(1);
      tstedge1   = Dr(1)*[-Inc*alp1(1)*fp1(1)+sum( an.*alp1.*fp1 );...
  	          -Inc*fm1(1) - sum( an.*fm1 )];
      tstedge2   = Dr(2)*[-sum( bn.*alp2.*fp2 ); sum( -bn.*fm2 )];
      tstedge    = [tstedge1,tstedge2]
   else
      disp('check frozen edge conditions');
      an         = -rn./Lam1;
      bn         = -tn./Lam2;
      Inc        = -1/Lam1(1);
      tstedge1   = [Dr(1)*(-Inc*alp1(1)*fp1(1)+sum(an.*alp1.*fp1));...
  	          -Dr(1)*( Inc*fm1(1) + sum(an.*fm1) );...
  	          Inc*alp1(1)-sum(an.*alp1);...
  	          Inc+sum(an)];
      tstedge2   = [-Dr(2)*sum( bn.*alp2.*fp2 );...
                    Dr(2)*sum( -bn.*fm2 );...
                   sum(bn.*alp2); sum(bn)];
      tstedge    = tstedge1-tstedge2
   end
 
   disp('test reverse wave');
   Rm  = -R'*T/T';
   Tm  = s*T;
   disp('expected R&T:')
   disp([Rm Tm]);
   %%
   if 1%%rerun program with thicknesses reversed to check
      [Rm2,Tm2] = SUB_RTstep_Galerkin(...
 			phys_vars,hh([2 1]),bc,NN,INC_SUB);
      disp('actual R&T:')
      disp([Rm2 Tm2]);
   end
end

function y = calc_res(Z2,Dr,hr,gamma)
%% y=calc_res(Z2,Dr,hr,gamma)=Res(1/f(K),gamma_n),
%% where gamma_n is a root of the dispersion relation
%% f=1/K/tanh(KH)-(Dr*K^4+lam-mr*mu);
%% Z2={lam,mu,H}.
lam   = Z2{1};
mu    = Z2{2};
H     = Z2{3};
mr    = hr;
%%
Gam   = Dr*gamma.^4+lam-mr*mu;
Gampr = Gam+4*Dr*gamma.^4;
denom = H*(Gam.^2.*gamma.^2-1)+Gampr;
y     = -gamma./denom;
