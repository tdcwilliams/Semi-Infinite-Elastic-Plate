function [Rp,Tp,Rm,Tm,y] = SUB_TMstep_Gal_comp(...
			phys_vars,hh,bc,NN,EE,SURGE)

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

INC_SUB  = 1;
do_test  = 0;
if nargin==0
   %% do some tests
   %% - print |R|&|T|, and check energy;
   do_test  = 1;
end

if ~exist('phys_vars')
   period      = 10;%% wave period [s]
   theta_inc   = 0;%% wave incident angle [degrees]
   H_dim       = 100;%% water depth [m]
   phys_vars   = {period,theta_inc,H_dim};
end
if ~exist('hh')
   %hh = [0 1];
   hh = [1 2];
   %hh = [2 1];
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

prams    = NDphyspram(0);%[E,g,rho_wtr,rho_ice,nu];
rho_wtr  = prams(3);
if ~exist('EE')
   EE = [prams(1),prams(1);
         prams(4),prams(4);
         prams(5),prams(5)];
else
   if prod(size(EE))==1
      EE = [EE,EE;
            prams(4),prams(4);
            prams(5),prams(5)];
   end
end

%%include surge or not
if ~exist('SURGE')
   SURGE = 1;
   %SURGE = 0;
end

%% tried to improve the number of roots needed
%% but not working at the moment - so don't recommend this option;

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
  [Rts,rts,HH,alpy,del0,L] = NDmakeZ2_sub(phys_vars,hh,Ninput,EE,rho_wtr);
else
  [Rts,rts,H,alpy,del0,L]  = NDmakeZ2_rel(phys_vars,hh,Ninput,EE,rho_wtr);
  HH                       = [H H];
end

%%want larger ice submergence on right:
DO_SWAP  = ( EE(2,1)*hh(1)>EE(2,2)*hh(2) );
if DO_SWAP
   hh    = fliplr(hh);
   EE    = fliplr(EE);
   Rts   = Rts([2 1]);
   rts   = rts([2 1]);
   HH    = HH([2 1]);
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
E1       = EE(1,1);
E2       = EE(1,2);
rho_ice1 = EE(2,1);
rho_ice2 = EE(2,2);
nu1      = EE(3,1);
nu2      = EE(3,2);
om       = 2*pi/phys_vars{1};

%%Compressional stuff for u problem:
mu1_lame = E1/2/(1+nu1);
mu2_lame = E2/2/(1+nu2);
Kc1_dim  = E1*hh(1)/(1-nu1^2);%%compressional rigidity ~ Pa*m ~ rho_wtr*om^2*L^2*L
Kc2_dim  = E2*hh(2)/(1-nu2^2);
Kc1      = Kc1_dim/(rho_wtr*om^2*L^2)/L;
Kc2      = Kc2_dim/(rho_wtr*om^2*L^2)/L;
%%
m1_dim   = rho_ice1*hh(1);
m2_dim   = rho_ice2*hh(2);
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
M1_rel   = rho_wtr*om^2*L^4;
   %% scale bending moment by this, since
   %% \int\sig_ij.(z-z_c).dz ~ Pa*m^2 ~ rho_wtr*om^2*L^2*L^2 = D1/L
M0_rel   = rho_wtr*om^2*L^3;
   %%scale zero-th moment by this, since
   %% \int\sig_ij.dz ~ Pa*m ~ rho_wtr*om^2*L^2*L = D1/L^2
J1_dim   = rho_ice1*hh(1)^3/12;
J2_dim   = rho_ice2*hh(2)^3/12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    tst_disprel   = ( Dr(j)*gg.^4+del0*[1;sigr(j)] ).*...
                     gg.*tanh(gg*HH(j))-1
  end
  return
end

%%
lam         = del0(1);
mu          = -del0(2);
nu_tilde1   = (1-nu1)*alpy^2;
nu_tilde2   = (1-nu2)*alpy^2;
H           = H2+sig2;
%%
gam0  = gam1(1);
alp0  = alp1(1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input1      = {gam1,alp1,H1;
               gam2,alp2,H2};
sigr        = [sig1/sig2,1];%H
nunu_tilde  = [nu_tilde1,nu_tilde2];
input2      = {del0,Dr,sigr,nunu_tilde};
%%
input3   = {kc1,Kc1;
            kc2,Kc2;
            zc1,zc2};
%%
M1    = 1;
M2    = 1;
NMM   = [Nterms,M1,M2];
DO_KC = 0;
%%
[MK,forcing,xtra,intrinsic_admittance] =...
   SUB_step_Gal_comp_kernel_forcing(input1,input2,input3,NMM,DO_KC);

%% solve integral eqn
uuB   = MK\forcing;

%xtra  = {rn_P1,tn_P2;
%         M_PM1,M_PM2;
%         M_u2r,M_u2t};
rn_unk   = xtra{1,1};
tn_unk   = xtra{1,2};
M_PM1B   = xtra{2,1}.';
M_PM2B   = xtra{2,2}.';
M_u2r    = xtra{3,1};
M_u2t    = xtra{3,2};
M_M0w    = xtra{4,1}.';
M_M1w    = xtra{4,2}.';
%%
j_inc    = 1:M1+M2+2;
j_unk    = M1+M2+2+(1:6);
jr_inc1  = 1:M1;
jc_inc1  = 1:M1;
jr_inc2  = 1:M2;
jc_inc2  = M1+(1:M2);
JC_comp  = M1+M2+(1:2);
%%
rr2                  = M_u2r*uuB;
tt2                  = M_u2t*uuB;
rr2(:,j_unk)         = rr2(:,j_unk)+rn_unk;
tt2(:,j_unk)         = tt2(:,j_unk)+tn_unk;
rr2(jr_inc1,jc_inc1) = rr2(jr_inc1,jc_inc1)+eye(M1);
tt2(jr_inc2,jc_inc2) = tt2(jr_inc2,jc_inc2)+eye(M2);
%%
if 0%%make testing easier
   jt    = [2 4 1 3];
   jt2   = 1:10;
   %Ainc_vec(jt)
   uuC   = [uuB(:,j_inc)*Ainc_vec(jt),uuB(:,j_unk)];
   rr2C  = [rr2(jt2,j_inc)*Ainc_vec(jt),rr2(jt2,j_unk)];
   tt2C  = [tt2(jt2,j_inc)*Ainc_vec(jt),tt2(jt2,j_unk)]
end

%%EDGE CONDITIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Dr(1)==0%%water on left:
   jd0      = [1 2 3 5];%%no more [S1 (1),Q1 (2),U1 (5)], & S2=0 (3)
                        %%keep Q2 (4) and U2 (6)
   Ninc     = M1+M2+1;%%no of incident waves
   jd1      = Ninc;%%no longer have incident compressive wave from left
   j_unk    = [Ninc+(1:2)];
   j_inc    = 1:Ninc;
   JC_comp  = Ninc;%%column for incident compressive wave from right
   j_dis2   = [jd1,M1+M2+2+jd0];
   %%
   rr2(:,j_dis2)  = [];
   tt2(:,j_dis2)  = [];
   
   %%bending moment eqn (RHS)
   M2B            = M_PM2B(:,2).'*tt2;
   M2B(jc_inc2)   = M2B(jc_inc2)+M_PM2B(jr_inc2,2).';%%'+' since M has even derivatives
   M2Bw           = M_M1w.'*rr2;
   M2Bw(jc_inc1)  = M2Bw(jc_inc1)+M_M1w(jr_inc1).';%%add pressure from incident f-g waves from left
   Medge          = M2B-M2Bw;

   if SURGE==1
      %%compression at edge (RHS): Bc_scat-Bc_inc = U2-2*Bc_inc = M0w/(1i*kc2*Kc2)
      M0w            = M_M0w.'*rr2;
      M0w(jc_inc1)   = M0w(jc_inc1)+ M_M0w(jr_inc1).';
      %%
      Medge(2,:)              = M0w/(1i*kc2*Kc2);
      Medge(2,JC_comp(end))   = Medge(2,JC_comp(end))+2;
      Medge(2,j_unk(end))     = Medge(2,j_unk(end))-1;%%u(0^+) column
   else
      %% set U2=0
      %% => the f-g problem is decoupled
      %% but the compressive problem depends on the f-g one:
      %% The compressive moment should still balance:
      %%  -2*Bc_inc = M0w/(1i*kc2*Kc2) (*)
      %%  => only 1 value of Bc_inc is prescribed from (*)
      Medge(2,j_unk(end))  = 1;
   end

   Q2B   = -Medge(:,j_unk)\Medge(:,j_inc);
   rn2   = rr2(:,j_inc);
   tn2   = tt2(:,j_inc);
   rn2   = rn2+rr2(:,j_unk)*Q2B;
   tn2   = tn2+tt2(:,j_unk)*Q2B;
   %%
   Bc_scat           = Q2B(2,:);%-Bc_inc;
   Bc_scat(JC_comp)  = Bc_scat(JC_comp)-1;
   if 0
      jt    = [2 4 3];
      jt2   = 1:10;
      tst_Q2B   = Q2B*Ainc_vec(jt);
      tst_rn2   = rn2(jt2,:)*Ainc_vec(jt);
      tst_tn2   = tn2(jt2,:)%*Ainc_vec(jt)
   end
   %%
   Rp = rn2(jr_inc1,jc_inc1)
   Tp = [Bc_scat(jc_inc1);
         tn2(jr_inc2,jc_inc1)]
   %%
   jc2   = [JC_comp,jc_inc2];
   Tm    = rn2(jr_inc1,jc2);
   Rm    = [Bc_scat(jc2);
            tn2(jr_inc2,jc2)];

   if SURGE==0
      Rm(:,1)  = [];%%inc compressive inc wave on RHS is not arbitrary
      Tm(:,1)  = [];
      Tp(1,:)  = [];%%=> compressive wave generated in RHS is not arbitrary
      Rm(1,:)  = [];
%     {'R',Rp(2,:)*[Ainc2]+Tm(2,:)*[Binc2];
%      'T',Tp(2,:)*[Ainc2]+Rm(2,:)*[Binc2]}
%  else
%     {'R',Rp*Ainc2+Tm(1,:)*[Bc_inc;Binc2];
%      'T',Tp(2,:)*Ainc2+Rm(2,:)*[Bc_inc;Binc2];
%      'Bc_scat',Tp(1,:)*Ainc2+Rm(1,:)*[Bc_inc;Binc2]}
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif bc==1%%free edges: col's of Medge corresp to [inc-waves (Ninc col's), psi1 psi2 U1 U2]
   jd0      = [1 3];%%S(0^-)=0 (col 1), S(0^+)=0 (col 3)
   Ninc     = M1+M2+2;%%no of incident waves
   j_unk    = [Ninc+(1:4)];
   j_inc    = 1:Ninc;
   JC_comp  = Ninc-1:Ninc;%%columns for incident compressive waves
   j_dis2   = Ninc+jd0;
   %%
   rr2(:,j_dis2)  = [];
   tt2(:,j_dis2)  = [];
   Medge          = zeros(4,size(rr2,2));
   
   %%Bending moment eqn (LHS)
   Medge(1,:)        = M_PM1B(:,2).'*rr2;%%bending moment: [L;R]
   Medge(1,jc_inc1)  = Medge(1,jc_inc1)+M_PM1B(jr_inc1,2).';
      %%need to add the LHS incident wave bending moment

   %%bending moment eqn (RHS)
   M2B            = M_PM2B(:,2).'*tt2;
   M2B(jc_inc2)   = M2B(jc_inc2)+M_PM2B(jr_inc2,2).';%%'+' since M has even derivatives
   M2Bw           = M_M1w.'*rr2;
   M2Bw(jc_inc1)  = M2Bw(jc_inc1)+M_M1w(jr_inc1).';%%add pressure from incident f-g waves from left
   Medge(2,:)     = M2B-M2Bw;

   if SURGE==1
      %%compression at edge (LHS):
      %% *U1=Ac_inc+Ac_scat=2*Ac_inc, or U1-2*Ac_inc=0
      Medge(3,j_unk(3))    = -1;
      Medge(3,JC_comp(1))  = 2;

      %%compression at edge (RHS): Bc_scat-Bc_inc = U2-2*Bc_inc = M0w/(1i*kc2*Kc2)
      M0w               = M_M0w.'*rr2;
      M0w(jc_inc1)      = M0w(jc_inc1)+ M_M0w(jr_inc1).';
      M0i               = 0*M0w;
      M0i(1,JC_comp(2)) = -2;
      M0i(1,j_unk(4))   = 1;
      Medge(4,:)        = M0w-1i*kc2*Kc2*M0i;
   else
      %%NB U1=0 => Ac_inc=0 so no compressive waves can exist in the LHS
      Medge(3,j_unk(3))    = 1;
      %% NB U2=0 => -2*Bc_inc = M0w/(1i*kc2*Kc2)
      %% => incident wave can only have one amplitude
      Medge(4,j_unk(4))    = 1;
   end
   %%
   Q2B   = -Medge(:,j_unk)\Medge(:,j_inc);
   rn2   = rr2(:,j_inc);
   tn2   = tt2(:,j_inc);
   rn2   = rn2+rr2(:,j_unk)*Q2B;
   tn2   = tn2+tt2(:,j_unk)*Q2B;
   %%
   Ac_scat              = Q2B(3,:);
   Ac_scat(JC_comp(1))  = Ac_scat(JC_comp(1))-1;
   Bc_scat              = Q2B(4,:);
   Bc_scat(JC_comp(2))  = Bc_scat(JC_comp(2))-1;
   if 0
      jt    = [2 4 3];
      jt2   = 1:10;
      tst_Q2B   = Q2B*Ainc_vec(jt);
      tst_rn2   = rn2(jt2,:)*Ainc_vec(jt);
      tst_tn2   = tn2(jt2,:)%*Ainc_vec(jt)
   end
   %%
   jc1   = [JC_comp(1),jc_inc1];
   Rp    = [Ac_scat(jc1);
            rn2(jr_inc1,jc1)]
   Tp    = [Bc_scat(jc1);
            tn2(jr_inc2,jc1)]
   %%
   jc2   = [JC_comp(2),jc_inc2];
   Tm    = [Ac_scat(jc2);
            rn2(jr_inc1,jc2)];
   Rm    = [Bc_scat(jc2);
            tn2(jr_inc2,jc2)];

   if SURGE==0
      Rp(:,1)  = [];%%no compressive inc wave can exist on LHS
      Tp(:,1)  = [];
      Rp(1,:)  = [];%%=> no compressive wave can be generated in LHS
      Tm(1,:)  = [];
      Rm(:,1)  = [];%%inc compressive inc wave on RHS is not arbitrary
      Tm(:,1)  = [];
      Tp(1,:)  = [];%%=> compressive wave generated in RHS is not arbitrary
      Rm(1,:)  = [];
%     {'R',Rp(2,:)*[Ainc2]+Tm(2,:)*[Binc2];
%      'T',Tp(2,:)*[Ainc2]+Rm(2,:)*[Binc2]}
%  else
%     {'R',Rp(2,:)*[Ac_inc;Ainc2]+Tm(2,:)*[Bc_inc;Binc2];
%      'T',Tp(2,:)*[Ac_inc;Ainc2]+Rm(2,:)*[Bc_inc;Binc2];
%      'Ac_scat',Rp(1,:)*[Ac_inc;Ainc2]+Tm(1,:)*[Bc_inc;Binc2];
%      'Bc_scat',Tp(1,:)*[Ac_inc;Ainc2]+Rm(1,:)*[Bc_inc;Binc2]}
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif bc==0%%frozen edges
   jd0      = [3 4];%%S(0^+)=S(0^-), \psi(0^+)=\psi(0^-), so remove "+"
   Ninc     = M1+M2+2;%%no of incident waves
   j_unk    = [Ninc+(1:4)];
   j_inc    = 1:Ninc;
   JC_comp  = Ninc-1:Ninc;%%columns for incident compressive waves
   j_dis2   = Ninc+jd0;
   %%
   rr2(:,j_dis2-2)  = rr2(:,j_dis2-2)+rr2(:,j_dis2);
   tt2(:,j_dis2-2)  = tt2(:,j_dis2-2)+tt2(:,j_dis2);
   rr2(:,j_dis2)  = [];
   tt2(:,j_dis2)  = [];

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Medge = zeros(4,size(rr2,2));
   if SURGE==1
      %% cty of horizontal displacement
      Medge(1,j_unk) = [0 dzc -1 1];
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% 0 moment (RHS): (1i*kc2*Kc2)*(Bc_scat-Bc_inc) = (1i*kc2*Kc2)*(U2-2*Bc_inc)
      M0i2              = 0*Medge(2,:);
      M0i2(JC_comp(2))  = -2;
      M0i2(j_unk(4))    = 1;  
      M0i2              = M0i2*(1i*kc2*Kc2);

      %% 0 moment (LHS): (1i*kc1*Kc1)*(Ac_inc-Ac_scat) = (1i*kc1*Kc1)*(2*Ac_inc-U1)
      M0i1              = 0*Medge(2,:);
      M0i1(JC_comp(1))  = 2;
      M0i1(j_unk(3))    = -1;  
      M0i1              = M0i1*(1i*kc1*Kc1);

      %% 0 moment from water:
      M0w            = M_M0w.'*rr2;
      M0w(jc_inc1)   = M0w(jc_inc1)+ M_M0w(jr_inc1).';

      %% full 0 moment eqn:
      Medge(2,:)  = M0w+M0i1-M0i2;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Bending moment (LHS)
      M1i1           = M_PM1B(:,2).'*rr2;%%bending moment: [L;R]
      M1i1(jc_inc1)  = M1i1(jc_inc1)+M_PM1B(jr_inc1,2).';
         %%need to add the LHS incident wave bending moment

      %% Bending moment (RHS)
      M1i2           = M_PM2B(:,2).'*tt2;
      M1i2(jc_inc2)  = M1i2(jc_inc2)+M_PM2B(jr_inc2,2).';
         %%need to add the RHS incident wave bending moment
         %%'+' since M has even derivatives

      %% Bending moment (water)
      M2Bw           = M_M1w.'*rr2;
      M2Bw(jc_inc1)  = M2Bw(jc_inc1)+M_M1w(jr_inc1).';
      Medge(3,:)     = M1i2-M1i1-M2Bw-dzc*M0i1;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   else
      %%SET U2=U2=0
      %%=>psi(0^-)=0
      %%get standing waves in both sides
      %%amplitudes from M0 & M1 equations
      Medge(1,j_unk(3)) = 1;
      Medge(2,j_unk(4)) = 1;
      Medge(3,j_unk(2)) = 1;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% w (LHS):
   M_W1           = M_PM1B(:,1).'*rr2;
   M_W1(jc_inc1)  = M_W1(jc_inc1)+M_PM1B(jr_inc1,1).';

   %% w (RHS):
   M_W2           = M_PM2B(:,1).'*tt2;
   M_W2(jc_inc2)  = M_W2(jc_inc2)+M_PM2B(jr_inc2,1).';

   %% w cts:
   Medge(4,:) = M_W2-M_W1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   Q2B   = -Medge(:,j_unk)\Medge(:,j_inc);
   rn2   = rr2(:,j_inc);
   tn2   = tt2(:,j_inc);
   rn2   = rn2+rr2(:,j_unk)*Q2B;
   tn2   = tn2+tt2(:,j_unk)*Q2B;
   %%
   Ac_scat              = Q2B(3,:);%%U(0^-)=Ac_inc+Ac_scat
   Ac_scat(JC_comp(1))  = Ac_scat(JC_comp(1))-1;
   Bc_scat              = Q2B(4,:);%%U(0^+)=Bc_inc+Bc_scat
   Bc_scat(JC_comp(2))  = Bc_scat(JC_comp(2))-1;
   %%
   jc1   = [JC_comp(1),jc_inc1];
   Rp    = [Ac_scat(jc1);
            rn2(jr_inc1,jc1)]
   Tp    = [Bc_scat(jc1);
            tn2(jr_inc2,jc1)]
   %%
   jc2   = [JC_comp(2),jc_inc2];
   Tm    = [Ac_scat(jc2);
            rn2(jr_inc1,jc2)];
   Rm    = [Bc_scat(jc2);
            tn2(jr_inc2,jc2)];

   if SURGE==0
      Rp(:,1)  = [];%%no compressive inc wave can exist on LHS
      Tp(:,1)  = [];
      Rp(1,:)  = [];%%=> no compressive wave can be generated in LHS
      Tm(1,:)  = [];
      Rm(:,1)  = [];%%inc compressive inc wave on RHS is not arbitrary
      Tm(:,1)  = [];
      Tp(1,:)  = [];%%=> compressive wave generated in RHS is not arbitrary
      Rm(1,:)  = [];
%     {'R',Rp(2,:)*[Ainc2]+Tm(2,:)*[Binc2];
%      'T',Tp(2,:)*[Ainc2]+Rm(2,:)*[Binc2]}
%  else
%     {'R',Rp(2,:)*[Ac_inc;Ainc2]+Tm(2,:)*[Bc_inc;Binc2];
%      'T',Tp(2,:)*[Ac_inc;Ainc2]+Rm(2,:)*[Bc_inc;Binc2];
%      'Ac_scat',Rp(1,:)*[Ac_inc;Ainc2]+Tm(1,:)*[Bc_inc;Binc2];
%      'Bc_scat',Tp(1,:)*[Ac_inc;Ainc2]+Rm(1,:)*[Bc_inc;Binc2]}
   end
end

%%REFLECTION AND TRANSMISSION COEFFICIENTS:
if DO_SWAP==1
   tmp = Rp;
   Rp  = Rm;
   Rm  = tmp;
   %%
   tmp = Tp;
   Tp  = Tm;
   Tm  = tmp;
   %%
   intrinsic_admittance = intrinsic_admittance([2 1 4 3]);
   %%
   y  = {Z_out,tn2,rn2,intrinsic_admittance};
else
   y  = {Z_out,rn2,tn2,intrinsic_admittance};
end

return;


if do_test==1%%test new results
   disp('R&T, |R|&|T|, |R|^2+s*|T|^2:');
   disp([R2,T2])
   disp([abs(R2),abs(T2)])
   disp(R2*R2'+s*T2*T2')
   %%
   fp1 = gam1.^2+nu_tilde1;
   fp2 = gam2.^2+nu_tilde2;

   %% New energy check:
   Ew_fac1  = -1/2/BG1(1);
   Ew_fac2  = -1/2/BG2(1);
   %Ew_fac1  = -alp1(1)/2/BG1(1);
   %Ew_fac2  = -alp2(1)/2/BG2(1);
   Ec_fac1  = kc1*Kc1;
   Ec_fac2  = kc2*Kc2;
   %%
   Ec_in    = Ec_fac1*abs(Ac_inc )^2 + Ec_fac2*abs(Bc_inc )^2%%compressive energy input
   Ec_out   = Ec_fac1*abs(Ac_scat)^2 + Ec_fac2*abs(Bc_scat)^2%%compressive energy output
   Ew_in    = Ew_fac1*abs(Ainc2  )^2 + Ew_fac2*abs(Binc2  )^2%%flex-grav energy input
   Ew_out   = Ew_fac1*abs(rn2(1) )^2 + Ew_fac2*abs(tn2(1) )^2%%flex-grav energy output
   tstE_new = [Ew_in+Ec_in,Ew_out+Ec_out]                    %%total energy: in=out
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
   fp1 = gam1.^2+nu_tilde1;
   fp2 = gam2.^2+nu_tilde2;
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
