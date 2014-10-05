function [R1,T1,R2,T2,Smat,y]=SUB_TMstep_Gal(...
			phys_vars,hh,bc,MM,NN,INC_SUB,EE,rho_wtr,DO_KC)
%% CALL: [R1,T1,R2,T2,Smat,y]=SUB_TMstep_Gal(...
%%			phys_vars,hh,bc,MM,NN,INC_SUB,EE,rho_wtr,DO_KC)
%% phys_vars = period or [period,theta_inc] or {period,theta_inc,H},
%%  where theta_inc is angle of incidence and z=-H is the sea floor
%%  (theta_inc=0, H=oo by default if they aren't entered);
%% hh = [h_left,h_right];
%% bc = 0 for frozen edges, bc = 1 for free edges;
%% MM = [M_left,M_right];
%% NN = [N_gegenbauers,N_eigenfxns];
%%
%% INC_SUB = 0 for draught = 0,
%%       or  1 for draught_j = rho_ice*h_j/rho_water; 
%%
%% EE = Young's modulus for ice if want to change from default
%%   or [YM_left,YM_right      - Young's moduli
%%       ID_left,ID_right      - ice densities
%%       PR_left,PR_right]     - Poisson's ratios
%% rho_wtr = water density
%% DO_KC = 1: try to improve convergence of kernel matrix
%%      or 0: leave as it was

do_test  = 0;%%do_test==1 => print |R|&|T|, and check energy
if nargin==0
   do_test  = 1;
   if 0
      hh = [0,1];
   elseif 0
      hh = [1,2];
      bc = 1;
   elseif 1
      hh = [1,2];
      bc = 0;
   end
end

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
if ~exist('DO_KC')
   DO_KC = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WANT ICE THAT IS SUBMERGED MORE DEEPLY ON THE RIGHT:
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%nu    = NDphyspram(5);
nunu        = EE(3,:);
nunu_tilde  = (1-nunu)*alpy^2;
%%
lam   = del0(1);
mu    = -del0(2);
sig2  = mu;
H     = H2+sig2;
gam0  = gam1(1);
alp0  = alp1(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%book-keeping
%%rows:
jr_inc1 = 1:M1;
jr_inc2 = 1:M2;

%%columns:
jc_inc1  = 1:M1;
jc_inc2  = M1+jr_inc2;
jinc     = 1:(M1+M2);
jc_P1    = M1+M2+(1:2);
jc_P2    = M1+M2+(3:4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET KERNEL MATRIX AND FORCING TERMS
input1   = {gam1,alp1,H1;
            gam2,alp2,H2};
input2   = {del0,Dr,sigr,nunu_tilde};
NMM      = [Npolys,M1,M2];
%%
[MK,forcing,xtra,intrinsic_admittance] =...
   SUB_step_Gal_kernel_forcing(input1,input2,NMM,INC_SUB,DO_KC);

%forcing  = {finc1,ME1;
%            finc2,ME2};
finc1 = forcing{1,1};
ME1   = forcing{1,2};
finc2 = forcing{2,1};
ME2   = forcing{2,2};

%xtra  = {rn_P1,tn_P2;
%         M_PM1,M_PM2;
%         M_u2r,M_u2t};
rn_P1 = xtra{1,1};
tn_P2 = xtra{1,2};
M_PM1 = xtra{2,1}.';
M_PM2 = xtra{2,2}.';
M_u2r = xtra{3,1};
M_u2t = xtra{3,2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_SWAP
   intrinsic_admittance = 1/intrinsic_admittance;
   y  = {{gam2/L,gam1/L},{alp2/L,alp1/L},alpy/L,[H2 H1]*L,intrinsic_admittance};
else
   y  = {{gam1/L,gam2/L},{alp1/L,alp2/L},alpy/L,[H1 H2]*L,intrinsic_admittance};
end

%%SOLVE INTEGRAL EQN:
uu =  MK\[finc1, finc2, ME1, ME2];
%%
rr                   = M_u2r*uu;
rr(jr_inc1,jc_inc1)  = rr(jr_inc1,jc_inc1)+eye(M1);
rr(:,jc_P1)          = rr(:,jc_P1)+rn_P1;%-...
%%
tt                   = M_u2t*uu;
tt(jr_inc2,jc_inc2)  = tt(jr_inc2,jc_inc2)+eye(M2);%Min1+Min2+(3:4)
tt(:,jc_P2)          = tt(:,jc_P2)+tn_P2;%...
%%

%%APPLY EDGE CONDITIONS:
if Dr(1)==0%%water on left:
  j_dis        = M1+M2+(1:3); %%remove these columns (w',S on left go; S=0 on right)
  j_unk        = M1+M2+1;     %%unknown constant (w' RHS) corresponds to this column
  rr(:,j_dis)  = [];
  tt(:,j_dis)  = [];
  %%
  Mom2            = M_PM2(:,2).'*tt;
  Mom2(1,jc_inc2) = Mom2(1,jc_inc2) + M_PM2(jr_inc2,2).';%%add bending moment from RHS incident waves
  Q2              = -Mom2(:,j_unk)\Mom2(:,jinc);
  %%
  rn  = rr(:,jinc)+rr(:,j_unk)*Q2;
  tn  = tt(:,jinc)+tt(:,j_unk)*Q2;
elseif bc==1%%free edges:
  j_dis        = M1+M2+[1 3]; %%remove these columns (S=0 on left and right)
  j_unk        = M1+M2+(1:2); %%unknown constants (w' on left and right) correspond to this column
  rr(:,j_dis)  = [];
  tt(:,j_dis)  = [];
  %%
  MM              = [M_PM1(:,2).'*rr; M_PM2(:,2).'*tt];
  MM(1,jc_inc1)   = MM(1,jc_inc1)+M_PM1(jr_inc1,2).'; %%add bending moment from LHS incident waves
  MM(2,jc_inc2)   = MM(2,jc_inc2)+M_PM2(jr_inc2,2).'; %%add bending moment from RHS incident waves
  %%
  QQ  = -MM(:,j_unk)\MM(:,jinc);
  rn  = rr(:,jinc)+rr(:,j_unk)*QQ;
  tn  = tt(:,jinc)+tt(:,j_unk)*QQ;
else%%frozen edges:
  j_dis        = M1+M2+(3:4); %%cts edges: w',S are the same on LHS and RHS so add forcings together
  j_unk        = M1+M2+(1:2); %%unknown constants (w',S on left) correspond to this column
  rr(:,j_unk)  = rr(:,j_unk)+rr(:,j_dis); %%add RHS to LHS
  tt(:,j_unk)  = tt(:,j_unk)+tt(:,j_dis); %%add RHS to LHS
  rr(:,j_dis)  = [];                      %%remove RHS
  tt(:,j_dis)  = [];                      %%remove RHS
  %%
  PM1             = M_PM1.'*rr;
  PM1(:,jc_inc1)  = PM1(:,jc_inc1)+M_PM1(jr_inc1,:).';%%add w,M from inc wave on LHS
  PM2             = M_PM2.'*tt;
  PM2(:,jc_inc2)  = PM2(:,jc_inc2)+M_PM2(jr_inc2,:).';%%add w,M from inc wave on RHS
  dPM             = PM2-PM1;
  %%
  Q1  = -dPM(:,j_unk)\dPM(:,jinc);
  rn  = rr(:,jinc)+rr(:,j_unk)*Q1;
  tn  = tt(:,jinc)+tt(:,j_unk)*Q1;
end

%%REFLECTION AND TRANSMISSION COEFFICIENTS:
jout1 = 1:M1;
jout2 = 1:M2;
R1    = rn(jr_inc1,jc_inc1);
T1    = tn(jr_inc2,jc_inc1);
R2    = tn(jr_inc2,jc_inc2);
T2    = rn(jr_inc1,jc_inc2);

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
   Rm = Smat(2,2);
   Tm = Smat(1,2);
   %%
   disp('R&T, |R|&|T|, |R|^2+s*|T|^2:');
   disp([Rp,Tp])
   disp([abs(Rp),abs(Tp)])
   %%
   s_ia  = y{end};
   disp(Rp*Rp'+s_ia*Tp*Tp')
   %E_tst = abs(Rp^2)+s_ia*abs(Tp^2)

   %%reverse waves:
   Rm2      = -Rp'*Tp/Tp';
   Tm2      = (1-abs(Rp^2))/Tp';
   rev_test = [Rm ,Tm ;
               Rm2,Tm2]
end
