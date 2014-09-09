function [R1,T1,R2,T2,Smat,y]=SUB_RTstep_Gal_manymodes(...
			phys_vars,hh,bc,MM,NN,INC_SUB,EE,rho_wtr)
%% CALL: [R1,T1,R2,T2,Smat,y]=SUB_RTstep_Gal_manymodes(...
%%			phys_vars,hh,bc,MM,NN,INC_SUB)
%% phys_vars = period or [period,theta_inc] or {period,theta_inc,H},
%%  where theta_inc is angle of incidence and z=-H is the sea floor
%%  (theta_inc=0, H=oo by default if they aren't entered);
%% hh = [h_left,h_right];
%% bc = 0 for frozen edges, bc = 1 for free edges;
%% MM = [M_left,M_right];
%% NN = [N_gegenbauers,N_eigenfxns];
%% INC_SUB = 0 for draught = 0,
%%  INC_SUB = 1 for draught_j = rho_ice*h_j/rho_water; 

do_test  = 0;%%do_test==1 => print |R|&|T|, and check energy

if ~exist('phys_vars');
   phys_vars   = {10,0,100};
end
if ~exist('hh');
   hh = [0 1];
end
if ~exist('bc');
   bc = 1;
end
if ~exist('MM')
   MM = [1 1];
end
if ~exist('NN')
   NN = [50 1000];
end
if ~exist('INC_SUB');
  INC_SUB   = 1;
end

prams    = NDphyspram(0);%[E,g,rho_wtr,rho_ice,nu];
if ~exist('EE')
   EE = [prams(1),prams(1);
         prams(4),prams(4);
         prams(5),prams(5)];
end
if ~exist('rho_wtr')
   rho_wtr  = prams(3);
end

%% DO NONDIMENSIONALIZATION:
Npolys   = NN(1);
if length(NN)==1
  decay_exp = 3;
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
  [Rts,rts,HH,alpy,del0,L] = NDmakeZ2_sub(phys_vars,hh,Ninput,EE,rho_wtr);
else
  [Rts,rts,H,alpy,del0,L]  = NDmakeZ2_rel(phys_vars,hh,Ninput);
  HH                       = [H H];
end

%%want ice that is submerged more deeply on the right:
sigsig   = EE(2,:).*hh/rho_wtr;
DO_SWAP  = ( sigsig(1)>sigsig(2) );
if DO_SWAP
  hh     = fliplr(hh);
  EE     = fliplr(EE);
  sigsig = fliplr(sigsig);
  Rts    = Rts([2 1]);
  rts    = rts([2 1]);
  HH     = HH([2 1]);
  MM     = fliplr(MM);
end

hr    = hh/hh(2);
Er    = EE(1,:)/EE(1,2);
sigr  = sigsig/sigsig(2);
Dr    = Er.*hr.^3;
M1    = MM(1);
M2    = MM(2);

%HH,20/L
gam1  = Rts{1};
alp1  = rts{1};
H1    = HH(1);
gam2  = Rts{2};
alp2  = rts{2};
H2    = HH(2);

%%
lam   = del0(1);
mu    = -del0(2);
nu    = NDphyspram(5);
nu1   = (1-nu)*alpy^2;
sig2  = mu;
H     = H2+sig2;
%%
gam0  = gam1(1);
alp0  = alp1(1);
%%
del1  = lam-sigr(1)*mu;
BGzz1 = calc_res({Dr(1),del1,H1},gam1).*gam1./alp1;
Lam1  = Dr(1)*gam1.^4+del1;
BGz1  = -Lam1.*BGzz1;
BG1   = -Lam1.*BGz1;
%%
del2  = lam-sigr(2)*mu;
BGzz2 = calc_res({Dr(2),del2,H2},gam2).*gam2./alp2;
Lam2  = Dr(2)*gam2.^4+del2;
BGz2  = -Lam2.*BGzz2;
BG2   = -Lam2.*BGz2;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%5th output:
intrinsic_admittance = ( BG1(1)/BG2(1) );
if DO_SWAP
   intrinsic_admittance = 1/intrinsic_admittance;
   y  = {{gam2/L,gam1/L},{alp2/L,alp1/L},alpy/L,[H2 H1]*L,intrinsic_admittance};
else
   y  = {{gam1/L,gam2/L},{alp1/L,alp2/L},alpy/L,[H1 H2]*L,intrinsic_admittance};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALC Fm0 & Fmr:
mvec  = (0:Npolys)';
alpC  = .5-1/3*INC_SUB*(sigr(1)~=1);%%=1/6 if submergence included;
			       %%=1/2 if no submergence.
kap1     = -1i*gam1*H2;
besJ1    = besselj(2*mvec'+alpC,kap1).';
c_left   = gamma(alpC)*(alpC+2*mvec).*(-1).^mvec;
c_rt1    = (2./kap1).^alpC./cosh(gam1*H1);
F1       = diag(c_left)*besJ1*diag(c_rt1);
%%
kap2  = -1i*gam2*H2;
besJ2 = besselj(2*mvec'+alpC,kap2).';
c_rt2 = (2./kap2).^alpC./cos(kap2);
F2    = diag(c_left)*besJ2*diag(c_rt2);
%%

%%MAIN KERNEL MATRIX:
MK = F2*diag(BG2)*F2.'+F1*diag(BG1)*F1.';

%%FORCING TERMS:
fm1   = gam1.^2-nu1;
E1    = [1+0*alp1,-Dr(1)*fm1];
fm2   = gam2.^2-nu1;
E2    = [1+0*alp2,-Dr(2)*fm2];

%%SOLVE INTEGRAL EQN:
jinc1 = 1:M1;
jinc2 = 1:M2;
jinc  = 1:(M1+M2);
finc1 = -F1(:,jinc1);
finc2 = F2(:,jinc2);
ME2   =  F2*diag(BGz2)*E2;
ME1   =  F1*diag(BGz1)*E1;
uu    =  MK\[finc1, finc2, ME1, ME2];
%%
rr                = 2*diag(BG1)*F1.'*uu;
rr(jinc1,jinc1)   = rr(jinc1,jinc1)+eye(M1);
rr(:,M1+M2+(1:2)) = rr(:,M1+M2+(1:2))-...
			2*diag(BGz1)*E1;%tstr=[rn,rr*[1;P1]]
%%
tt                   = -2*diag(BG2)*F2.'*uu;
tt(jinc2,M1+jinc2)   = tt(jinc2,M1+jinc2)+eye(M2);%Min1+Min2+(3:4)
tt(:,M1+M2+(3:4))    = tt(:,M1+M2+(3:4))+...
			2*diag(BGz2)*E2;%tstt=[tn tt*[1;P1]]
%%
M_PM1 = E1.'*diag(-1./Lam1);
M_PM2 = E2.'*diag(-1./Lam2);

%%APPLY EDGE CONDITIONS:
if Dr(1)==0%%water on left:
  j_dis        = M1+M2+(1:3);
  j_mat        = M1+M2+1;
  rr(:,j_dis)  = [];
  tt(:,j_dis)  = [];
  %%
  Mom2               = M_PM2(2,:)*tt;
  Mom2(1,M1+jinc2)   = Mom2(1,M1+jinc2) + M_PM2(2,jinc2);
  Q2                 = -Mom2(j_mat)\Mom2(jinc);%M2(1)/M2(2);
  %%
  rn  = rr(:,jinc)+rr(:,j_mat)*Q2;
  tn  = tt(:,jinc)+tt(:,j_mat)*Q2;
elseif bc==1%%free edges:
  j_dis        = M1+M2+[1 3];
  j_mat        = M1+M2+(1:2);
  rr(:,j_dis)  = [];
  tt(:,j_dis)  = [];
  %%
  MM              = [M_PM1(2,:)*rr; M_PM2(2,:)*tt];
  MM(1,jinc1)     = MM(1,jinc1)+M_PM1(2,jinc1);
  MM(2,M1+jinc2)  = MM(2,M1+jinc2)+M_PM2(2,jinc2);
  %%
  QQ  = -MM(:,j_mat)\MM(:,jinc);
  rn  = rr(:,jinc)+rr(:,j_mat)*QQ;
  tn  = tt(:,jinc)+tt(:,j_mat)*QQ;
else%%frozen edges:
  j_dis        = M1+M2+(3:4);
  j_mat        = M1+M2+(1:2);
  rr(:,j_mat)  = rr(:,j_mat)+rr(:,j_dis);
  tt(:,j_mat)  = tt(:,j_mat)+tt(:,j_dis);
  rr(:,j_dis)  = [];
  tt(:,j_dis)  = [];
%    rr(:,2:3)=rr(:,2:3)+rr(:,4:5);
%    rr(:,4:5)=[];
%    tt(:,2:3)=tt(:,2:3)+tt(:,4:5);
%    tt(:,4:5)=[];
  %%
  PM1             = M_PM1*rr;
  PM1(:,jinc1)    = PM1(:,jinc1)+M_PM1(:,jinc1);
  PM2             = M_PM2*tt;
  PM2(:,M1+jinc2) = PM2(:,M1+jinc2)+M_PM2(:,jinc2);
  dPM             = PM2-PM1;
  %%
  Q1  = -dPM(:,j_mat)\dPM(:,jinc);
  rn  = rr(:,jinc)+rr(:,j_mat)*Q1;
  tn  = tt(:,jinc)+tt(:,j_mat)*Q1;
end

%%REFLECTION AND TRANSMISSION COEFFICIENTS:
jout1 = 1:M1;
jout2 = 1:M2;
R1    = rn(jout1,jinc1);
T1    = tn(jout2,jinc1);
R2    = tn(jout2,M1+jinc2);
T2    = rn(jout1,M1+jinc2);

if DO_SWAP
  tmp = R2;
  R2  = R1;
  R1  = tmp;
  %%
  tmp = T2;
  T2  = T1;
  T1  = tmp;
end
Smat  = [R1(1) T2(1);T1(1) R2(1)];

if do_test==1
   Rp = Smat(1,1);
   Tp = Smat(2,1);
   Rm = Smat(1,2);
   Tm = Smat(2,2);
   %%
   s_ia  = y{end};
   E_tst = abs(Rp^2)+s_ia*abs(Tp^2)

   %%reverse waves:
   Rm2      = -Rp'*Tp/Tp';
   Tm2      = (1-abs(Rp^2))/Tp';
   rev_test = [Rm ,Tm ;
               Rm2,Tm2]
end

function y=calc_res(Z2,gamma)
%% y=calc_res(Z2,hr,gamma)=Res(1/f(K),gamma_n),
%% where gamma_n is a root of the dispersion relation
%% f=1/K/tanh(KH)-(Dr*K^4+del);
%% Z2={Dr,lam-mr*mu,H}.
Dr    = Z2{1};
del   = Z2{2};
H     = Z2{3};
%%
Gam   = Dr*gamma.^4+del;%%now include mr in mu
Gampr = Gam+4*Dr*gamma.^4;
denom = H*(Gam.^2.*gamma.^2-1)+Gampr;
y     = -gamma./denom;
