function [R,T,y] = SUB_RTstep_Gal_comp_v0(...
			phys_vars,hh,bc,NN,INC_SUB,do_test,Mlog)
%% same as Williams & Porter (2009),
%% but changed to make notation closer to Williams (2014) 
%%
%% calc's scattering coeffiecients for an ice step
%% CALL: [R,T,y]=...
%%   SUB_RTstep_Gal_comp_v0(phys_vars,hh,bc,NN,INC_SUB)
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
end
if ~exist('bc')
   bc = 1;%% free edge conditions;
   %bc = 0;%% frozen edge conditions;
end
if ~exist('NN')
   NN = [10 1000];%% [N_poly, N_roots];
end
if ~exist('INC_SUB')
   INC_SUB = 1;%% include submergence or not;
end

if nargin==0
   %% do some tests
   %% - print |R|&|T|, and check energy;
   do_test  = 1;
else
   if ~exist('do_test')
      do_test  = 0;
   end
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
  [Rts,rts,HH,alpy,del0,L] = NDmakeZ2_sub(phys_vars,hh,Ninput);
else
  [Rts,rts,H,alpy,del0,L]  = NDmakeZ2_rel(phys_vars,hh,Ninput);
  HH                       = [H H];
end

%%want larger ice thickness on right:
DO_SWAP  = ( hh(1)>hh(2) );
if DO_SWAP==0
   Ainc2 = 1;%%wave from left
   Binc2 = 0;
else
   hh    = fliplr(hh);
   Rts   = Rts([2 1]);
   rts   = rts([2 1]);
   HH    = HH([2 1]);
   %%
   Ainc2 = 0;
   Binc2 = 1;%%wave from right now
end

if ~DO_SWAP
   Z_out = {{Rts{1}/L,Rts{2}/L},{rts{1}/L,rts{2}/L},...
             HH*L,alpy/L,del0,L};
else
   Z_out = {{Rts{2}/L,Rts{1}/L},{rts{2}/L,rts{1}/L},...
             HH([2 1])*L,alpy/L,del0,L};
end


%%Elastic constants:
nu       = NDphyspram(5);
E1       = NDphyspram(1);%%E the same on both sides
E2       = NDphyspram(1);
mu1_lame = E1/2/(1+nu);
mu2_lame = E2/2/(1+nu);

%%
hr = hh/hh(2);
Dr = [E1/E2,1].*hr.^3;
%%
gam1  = Rts{1};%gam1(1:10)
alp1  = rts{1};
H1    = HH(1);
gam2  = Rts{2};
alp2  = rts{2};
H2    = HH(2);

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
BGzz1 = calc_res({lam,mu,H1},hr(1),gam1).*gam1./alp1;%%BG_0^(0)
Lam1  = Dr(1)*gam1.^4+lam-hr(1)*mu;
BGz1  = -Lam1.*BGzz1;%%BG_0^(1); why '-'? (in W&P paper, but seems wrong)
BGz1B = Lam1.*BGzz1;%%fix minus sign in BG_0^(1)
BG1   = -Lam1.*BGz1;%%BG_0^(2)
%%
BGzz2 = calc_res({lam,mu,H2},hr(2),gam2).*gam2./alp2;
Lam2  = Dr(2)*gam2.^4+lam-hr(2)*mu;
BGz2  = -Lam2.*BGzz2;%%why '-'? (in W&P paper, but seems wrong)
BGz2B = Lam2.*BGzz2;%%fix minus sign
BG2   = -Lam2.*BGz2;
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
kap2  = -i*gam2*H2;
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
  
%%FORCING TERMS:
fm1   =  gam1.^2-nu1;
fm2   =  gam2.^2-nu1;
E1    =  [1+0*alp1,-Dr(1)*fm1];
E2    =  [1+0*alp2,-Dr(2)*fm2];
E1B   =  [-1+0*alp1,Dr(1)*fm1];%% S=-L^2*D/D1*w_xxx, psi=-w_x
E2B   =  [-1+0*alp2,Dr(2)*fm2];

%%SOLVE INTEGRAL EQN:
%Ainc  = -1i;
Ainc  = 1;
ME2   = F2*diag(1i*BGz2)*E2;
ME1   = F1*diag(1i*BGz1)*E1;
uu    = MK\[Ainc*F1(:,1),ME1,ME2];%%(38a)
ME2B  = F2*diag(1i*BGz2B)*E2B;
ME1B  = F1*diag(1i*BGz1B)*E1B;
uuB   = MK\[Ainc2*F1(:,1)-Binc2*F2(:,1),ME1B,ME2B];
%%
rr          = 2*diag(BG1)*F1.'*uu;%%in (24) of W&P: rr=1i*R_0(\alpha)
rr(1)       = rr(1)+1i*Ainc;
rr(:,2:3)   = rr(:,2:3)-2*diag(BGz1)*E1;
%%
rr2          = -2i*diag(BG1)*F1.'*uuB;%%now same as (24) of W&P (but Q vector -> -Q)
rr2(1)       = rr2(1)+Ainc2;%%now same as (24) of W&P
rr2(:,2:3)   = rr2(:,2:3)+2i*diag(BGz1B)*E1B;
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
M_PM1B   = diag([-1,1])*E1B.'*diag(1./Lam1);%% M=-L*Dr*w_xx
M_PM2B   = diag([-1,1])*E2B.'*diag(1./Lam2);%%

%%APPLY EDGE CONDITIONS:
COMP_OLD = 0;
if Dr(1)==0%%water on left:
   j_dis        = 2:4;
   rr(:,j_dis)  = [];
   tt(:,j_dis)  = [];
   rr2(:,j_dis) = [];
   tt2(:,j_dis) = [];
   %%
   M2  = M_PM2(2,:)*tt;
   Q2  = -M2(1)/M2(2);
   M2B      = M_PM2B(2,:)*tt2;
   M2B(1)   = M2B(1)+Binc2*M_PM2B(2);%%'+' since M has even derivatives
   Q2B      = -M2B(1)/M2B(2);
   %%
   rn  = rr*[1;Q2];
   tn  = tt*[1;Q2];
   rn2 = rr2*[1;Q2B];
   tn2 = tt2*[1;Q2B];
   %%
   v_unk       = [1;Q2];
   v_unk2      = eye(3,1);
   v_unk2(3)   = Q2B;
   uu(:,j_dis) = [];
elseif bc==1%%free edges:
   rr(:,[2 4]) = [];
   tt(:,[2 4]) = [];%% S(L)=S(R)=0
   rr2(:,[2 4])   = [];
   tt2(:,[2 4])   = [];%% S(L)=S(R)=0
   %%
   MM    = [M_PM1(2,:)*rr; M_PM2(2,:)*tt];
   MM(1) = MM(1)+1i*Ainc*M_PM1(2,1);
   QQ    = -MM(:,2:3)\MM(:,1);%%[w_x L,w_x R]
   %%
   MM2      = [M_PM1B(2,:)*rr2; M_PM2B(2,:)*tt2];%%bending moment: [L;R]
   MM2(1,1) = MM2(1,1)+Ainc2*M_PM1B(2,1);%%need to add the LHS incident wave bending moment
   MM2(2,1) = MM2(2,1)+Binc2*M_PM2B(2,1);%%need to add the RHS incident wave bending moment
   QQ2      = -MM2(:,2:3)\MM2(:,1);%%[w_x L,w_x R]
   %%
   v_unk          = eye(5,1);
   v_unk([3 5])   = QQ;
   rn             = rr*[1;QQ];
   tn             = tt*[1;QQ];
   %%
   v_unk2         = eye(5,1);
   v_unk2([3 5])  = QQ2;
   rn2            = rr2*[1;QQ2];
   tn2            = tt2*[1;QQ2];
   if 0
      %[M_PM1(2,:)*rr; M_PM2(2,:)*tt]/1i
      %[M_PM1(2,:)*rr2; M_PM2(2,:)*tt2]
      MM/1i,MM2
      QQ,QQ2,return
   end
else%%frozen edges:
   rr(:,2:3)   = rr(:,2:3)+rr(:,4:5);%% [-S,w_x](L)=[-S,w_x](R) so eliminate (R)
   rr(:,4:5)   = [];
   rr2(:,2:3)  = rr2(:,2:3)+rr2(:,4:5);%% [-S,w_x](L)=[-S,w_x](R) so eliminate (R)
   rr2(:,4:5)  = [];
   %%
   tt(:,2:3)   = tt(:,2:3)+tt(:,4:5);
   tt(:,4:5)   = [];
   tt2(:,2:3)  = tt2(:,2:3)+tt2(:,4:5);
   tt2(:,4:5)  = [];
   %%
   PM1      = M_PM1*rr;
   PM1(:,1) = PM1(:,1)+1i*Ainc*M_PM1(:,1);
   PM2      = M_PM2*tt;
   dPM      = PM2-PM1;
   Q1       = -dPM(:,2:3)\dPM(:,1);
   %%
   PM1B        = M_PM1B*rr2;
   PM1B(:,1)   = PM1B(:,1)+Ainc2*M_PM1B(:,1);
   PM2B        = M_PM2B*tt2;
   PM2B(:,1)   = PM2B(:,1)+Binc2*M_PM2B(:,1);
   dPMB        = PM2B-PM1B;
   Q1B         = -dPMB(:,2:3)\dPMB(:,1);
   %%
   v_unk = [1;Q1;Q1];
   rn    = rr*[1;Q1];
   tn    = tt*[1;Q1];
   %%
   v_unk2   = [1;Q1B;Q1B];
   rn2      = rr2*[1;Q1B];
   %rn2      = rr2*[0;Q1B];'hey: test!'
   tn2      = tt2*[1;Q1B];
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
   s2 = ( BG1(1)/BG2(1) );%intrinsic admittance
   y2 = {Z_out,rn2,tn2,s2};
else
   R2 = tn2(1)/Binc2;
   T2 = rn2(1)/Binc2;
   s2 = ( BG2(1)/BG1(1) );%intrinsic admittance
   y2 = {Z_out,tn2,rn2,s2};
   %%
   %v_unk2   = v_unk2([1 4 5 2 3]);
end
R  = -1i*rn(1)/Ainc;
T  = -1i*tn(1)/Ainc;
s  = ( BG1(1)/BG2(1) );%intrinsic admittance
y  = {Z_out,rn,tn,s};

if COMP_OLD
   if DO_SWAP==0%%test R,T for LHS wave
      DO_SWAP
      [R,R2]
      [T,T2]
      return
   else
      DO_SWAP
      Rm  = -R'*T/T';
      Tm  = (1-abs(R)^2)/T';
      [Rm,R2]
      [Tm,T2]
      return
   end
end

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
   disp(R2*R2'+s2*T2*T2')
   %%
   fp1 = gam1.^2+nu1;
   fp2 = gam2.^2+nu1;


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

   Vunk_chk = [S_lhs;psi_lhs;S_rhs;psi_rhs];

   if hh(1)==0
      disp('test definitions of unknown constants')
      [v_unk2(2:3),Vunk_chk(3:4)]
      disp('check free edge conditions');
      tstedge    = wpms2(3:4)
   elseif bc==1
      disp('test definitions of unknown constants')
      [v_unk2(2:5),Vunk_chk]
      disp('check free edge conditions');
      Inc        = -1/Lam1(1);
      tstedge    = [wpms1(3:4),wpms2(3:4)]
   else
      disp('test definitions of unknown constants')
      [v_unk2(2:5),Vunk_chk]
      disp('check frozen edge conditions');
      tstedge    = [wpms1,wpms2]
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

function y = calc_res(Z2,hr,gamma)
%% y=calc_res(Z2,hr,gamma)=Res(1/f(K),gamma_n),
%% where gamma_n is a root of the dispersion relation
%% f=1/K/tanh(KH)-(Dr*K^4+lam-mr*mu);
%% Z2={lam,mu,H}.
lam   = Z2{1};
mu    = Z2{2};
H     = Z2{3};
Dr    = hr^3;
mr    = hr;
%%
Gam   = Dr*gamma.^4+lam-mr*mu;
Gampr = Gam+4*Dr*gamma.^4;
denom = H*(Gam.^2.*gamma.^2-1)+Gampr;
y     = -gamma./denom;
