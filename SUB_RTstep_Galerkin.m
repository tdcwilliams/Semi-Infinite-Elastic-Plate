function [R,T,y] = SUB_RTstep_Galerkin(...
			phys_vars,hh,bc,NN,INC_SUB,youngs,rho_wtr,do_test,Mlog)
%% calc's scattering coeffiecients for an ice step
%% CALL: [R,T,y]=...
%%   SUB_RTstep_Galerkin(phys_vars,hh,bc,NN,INC_SUB)
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
   %hh = [1 2];
   %hh = [2 1];
   hh = [0 1];
   %hh = [1 0];
end
if ~exist('bc')
   bc = 1;%% free edge conditions;
   %bc = 0;%% frozen edge conditions;
end
if ~exist('NN')
   NN = [50 1000];%% [N_poly, N_roots];
end
if ~exist('INC_SUB')
   INC_SUB = 1;%% include submergence or not;
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

do_test  = 0;
if nargin==0
   %% do some tests
   %% - print |R|&|T|, and check energy;
   do_test  = 1;
end

SEP_FXN  = 0;%%use separate function for kernel matrix and focing vectors;

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
   [Rts,rts,HH,alpy,del0,L] = NDmakeZ2_sub(phys_vars,hh,Ninput,EE,rho_wtr);
else
   [Rts,rts,H,alpy,del0,L]  = NDmakeZ2_rel(phys_vars,hh,Ninput,EE,rho_wtr);
   HH                       = [H H];
end

%%want larger ice thickness on right:
sigsig   = EE(2,:).*hh/rho_wtr;
DO_SWAP  = ( sigsig(1)>sigsig(2) );
%DO_SWAP  = ( hh(1)>hh(2) );
if DO_SWAP
  hh     = fliplr(hh);
  %youngs = fliplr(youngs);
  EE     = fliplr(EE);
  sigsig = fliplr(sigsig);
  Rts    = Rts([2 1]);
  rts    = rts([2 1]);
  HH     = HH([2 1]);
end

if ~DO_SWAP
   Z_out = {{Rts{1}/L,Rts{2}/L},{rts{1}/L,rts{2}/L},...
             HH*L,alpy/L,del0,L};
else
   Z_out = {{Rts{2}/L,Rts{1}/L},{rts{2}/L,rts{1}/L},...
             HH([2 1])*L,alpy/L,del0,L};
end
%Z_out{1}{1}(1:10),Z_out{1}{2}(1:10)

hr    = hh/hh(2);
Er    = EE(1,:)/EE(1,2);
sigr  = sigsig/sigsig(2);
Dr    = Er.*hr.^3;
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
      gg          = Rts{j};
      tst_disprel = ( Dr(j)*gg.^4+del0*[1;sigr(j)] ).*...
                      gg.*tanh(gg*HH(j))-1
   end
   return
end

%%
lam         = del0(1);
mu          = -del0(2);
nunu        = EE(3,:);%NDphyspram(5);
nunu_tilde  = (1-nunu)*alpy^2;
sig2        = mu;
H           = H2+sig2;
%%
gam0  = gam1(1);
alp0  = alp1(1);
%%

if SEP_FXN==0
   %%Coefficients of Green's fxn expansion
   %BGzz1 = calc_res({lam,mu,H1},hr(1),gam1).*gam1./alp1;
   %Lam1  = Dr(1)*gam1.^4+lam-hr(1)*mu;
   del1  = lam-sigr(1)*mu;
   BGzz1 = calc_res({Dr(1),del1,H1},gam1).*gam1./alp1;
   Lam1  = Dr(1)*gam1.^4+del1;
   BGz1  = -Lam1.*BGzz1;
   BG1   = -Lam1.*BGz1;
   %%
   %BGzz2 = calc_res({lam,mu,H2},hr(2),gam2).*gam2./alp2;
   %Lam2  = Dr(2)*gam2.^4+lam-hr(2)*mu;
   del2  = lam-sigr(2)*mu;
   BGzz2 = calc_res({Dr(2),del2,H2},gam2).*gam2./alp2;
   Lam2  = Dr(2)*gam2.^4+del2;
   BGz2  = -Lam2.*BGzz2;
   BG2   = -Lam2.*BGz2;
   %%
   %tstBG1   = {BG1(1),H1,lam, mu,hr(1),Dr(1),gam1(1),L}
   %tstBG2   = {BG2(1),H2,lam, mu,hr(2),Dr(2),gam2(1),L}
   %%
   if 0
      phys_vars
      hh
      NN
      Rts{1}(1:5)
      Rts{2}(1:5)
      BG1(1:5)
      BG2(1:5)
      pause
   end

   %% CALC Fm0 & Fmr:
   mvec  = (0:Nterms)';
   alpC  = .5-1/3*INC_SUB*(hr(1)~=1);%%=1/6 if submergence included;
                                  %%=1/2 if no submergence.
   if 0%%old school
      kap1  = -1i*gam1*H2;
      besJ1 = besselj(2*mvec'+alpC,kap1).';
      kap2  = -1i*gam2*H2;
      besJ2 = besselj(2*mvec'+alpC,kap2).';
   else%%new school
      kap1        = -1i*gam1*H2;
      Nu1         = 2*mvec'+alpC;
      [NU1,KAP1]  = meshgrid(Nu1,kap1);
      besJ1       = besselj(NU1,KAP1).';
      %%
      kap2        = -1i*gam2*H2;
      [NU2,KAP2]  = meshgrid(Nu1,kap2);
      besJ2       = besselj(NU2,KAP2).';
   end
   c_left   = gamma(alpC)*(alpC+2*mvec).*(-1).^mvec;
   c_rt1    = (2./kap1).^alpC./cosh(gam1*H1);
   F1       = diag(c_left)*besJ1*diag(c_rt1);
   %%
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
   if ~FAST_KERNEL%%basic way:
      MK  = F2*diag(i*BG2)*F2.'+F1*diag(i*BG1)*F1.';
   elseif 1%%use correction:
      MK     = F2*diag(i*BG2)*F2.'+F1*diag(i*BG1)*F1.';
      nnvec  = (1:Nroots)';
      BGap   = 1i/pi./nnvec;
      %%
      gam1ap = 1i*nnvec/H1;
      kap1   = -i*H2*gam1ap;
      besJ1  = besselj(2*mvec'+alpC,kap1).';
      c_left = gamma(alpC)*(alpC+2*mvec).*(-1).^mvec;
      c_rt1  = (2./kap1).^alpC./cosh(gam1ap*H1);
      F1ap   = diag(c_left)*besJ1*diag(c_rt1);
      %%
      gam2ap = 1i*nnvec/H2;
      kap2   = -1i*H2*gam2ap;
      besJ2  = besselj(2*mvec'+alpC,kap2).';
      c_left = gamma(alpC)*(alpC+2*mvec).*(-1).^mvec;
      c_rt2  = (2./kap2).^alpC./cosh(gam2ap*H2);
      F2ap   = diag(c_left)*besJ2*diag(c_rt2);
      %%
      MK  = MK-F2ap*diag(1i*BGap)*F2ap.'+F1ap*diag(1i*BGap)*F1ap.';
      %%
      Nint      = 200;
      [tau,wgJ] = OP_numint_jacobi(alpC-.5,0,Nint);
      tt        = (1+tau)/2;
      ww        = wgJ.*( (3+tau)/4 ).^(alpC-.5);
      Npolys    = 2*Nterms;
      IP        = OP_inprod_gegenbauer(tt,ww,alpC,Npolys);
      IP        = IP(1:2:Npolys+1,:);
      %%
      zz           = -H+H2*tt;
      [Z0,Z]       = meshgrid(zz,zz);
      [T0,T]       = meshgrid(tt,tt);
      Jp           = find(abs(Z-Z0)>1e-9);
      Gcorr1       = diag(log(pi/H1+0*tt)/2/pi);
      Gcorr1(Jp)   = 1/2/pi*log(abs(1-exp(i*pi*(Z(Jp)-Z0(Jp))/H1))./abs(T(Jp)-T0(Jp)));
      Gcorr1       = Gcorr1+1/2/pi*log(abs(1-exp(1i*pi*(Z+Z0+2*H)/H1))./abs(T+T0));
      %%
      Gcorr2       = diag(log(pi/H2+0*tt)/2/pi);
      Gcorr2(Jp)   = 1/2/pi*log(abs(1-exp(i*pi*(Z(Jp)-Z0(Jp))/H2))./abs(T(Jp)-T0(Jp)));
      Gcorr2       = Gcorr2+1/2/pi*log(abs(1-exp(i*pi*(Z+Z0+2*H)/H2))./abs(T+T0));
      %[size(Gcorr1),size(Gcorr2),size(IP),size(MK)]
      MK  = MK+Mlog(mvec+1,mvec+1)/pi+IP*(Gcorr1+Gcorr2)*IP';
   else%%do more accurately, using recurrence relation:
      %MK=F2*diag(BG2)*F2.'+F1*diag(BG1)*F1.';
      tst00  = [];%MK(1,1);
      %%
      MK        = zeros(Nterms+1,Nterms+1);
      MK(1,1)   = F2(1,:)*diag(BG2)*F2(1,:).'+...
                    F1(1,:)*diag(BG1)*F1(1,:).';
      %tst00=[tst00,MK(1,1)];

      %%IMPROVE ACCURACY OF THIS TERM:
      fac0   = ( 2^alpC*alpC*gamma(alpC) )^2;
      fac1   = pi*H2/H1;
      Nlast  = 1e6;
      N10    = 1+2*(hh(1)>0);
      N1     = length(gam1)-N10;
      if 0%%check asympotics of the sums:
         nn   = (1:N1)';
         jt   = N1-10:N1;
         kap  = nn*pi*H2/H1;
         tst  = BG1(N10+nn).*(F1(1,N10+nn).').^2;
         tst1 = [tst,...
                1i/(H1/H2)*fac0/pi.*(1+cos(2*kap-alpC*pi-pi/2))./...
                   kap.^(2*alpC+2),...
                1i*fac0*fac1/pi^2*(1+cos(2*kap-alpC*pi-pi/2))./...
                   kap.^(2+2*alpC)];tst1(jt,:)
      end
      kap_xtra  = fac1*(N1+1:Nlast)';
      S         = sum((1+cos(2*kap_xtra-alpC*pi-pi/2))./kap_xtra.^(2+2*alpC));
      MK(1,1)   = MK(1,1)+i*fac0*fac1/pi^2*S;
      %%
      fac2   = pi;
      N20    = 1+2*(hh(2)>0);
      N2     = length(gam2)-N20;
      if 0%%check asympotics of the sums:
         nn   = (1:N2)';
         kap  = nn*pi*H2/H2;
         tst  = BG2(N20+nn).*(F2(1,N20+nn).').^2;
         jt   = N2-10:N2;
         tst2 = [tst,...
                1i*fac0*fac2/pi^2*(1+cos(2*kap-alpC*pi-pi/2))./...
                   kap.^(2+2*alpC)];tst2(jt,:)
         return;
      end
      kap_xtra  = fac2*(N2+1:Nlast)';
      S         = sum((1+cos(2*kap_xtra-alpC*pi-pi/2))./kap_xtra.^(2+alpC));
      MK(1,1)   = MK(1,1)+i*fac0*fac2/pi^2*S;
      %tst00=[tst00,MK(1,1)],return;
      clear kap_xtra;

      %%NOW GET REST OF MATRIX USING R.R.:
      if Nterms>0
         jj0        = (1:Nterms)';
         besJ1bar   = besJ1(jj0,:)+besJ1(jj0+1,:);
         %% =2/kap1*(2*m+alpC-1)*J_{2*m+alpC-1}(kap1)
         %% NB 1<=m<=Nterms
         F1bar      = diag(c_left(1+jj0))*besJ1bar*diag(c_rt1.*BG1);
         %%
         besJ2bar   = besJ2(jj0,:)+besJ2(jj0+1,:);
         %% =2/kap2*(2*m+alpC-1)*J_{2*m+alpC-1}(kap2)
         F2bar   = diag(c_left(1+jj0))*besJ2bar*diag(c_rt2.*BG2);
         MKbar   = F1*F1bar.'+F2*F2bar.';%%Nterms+1 x Nterms matrix
         %%
         for j=2:Nterms+1
            MK(j,1)  = MKbar(1,j-1)-c_left(j)/c_left(j-1)*MK(j-1,1);
         end
         for j=2:Nterms+1
            MK(:,j)  = MKbar(:,j-1)-c_left(j)/c_left(j-1)*MK(:,j-1);
         end
      end
   end

   %%FORCING TERMS:
   fm1   =  gam1.^2-nunu_tilde(1);
   E1    =  [1+0*alp1,-Dr(1)*fm1];
   fm2   =  gam2.^2-nunu_tilde(2);
   E2    =  [1+0*alp2,-Dr(2)*fm2];

   %ME2   = -F2*diag(1i*BGz2)*E2;
   %ME1   = -F1*diag(1i*BGz1)*E1;
   ME2   = F2*diag(BGz2)*E2;
   ME1   = F1*diag(BGz1)*E1;

   %%forcing from incident wave;
   finc1 = -F1(:,1);

   %%scattering caused by "forced oscillations" (\bfP_0 & \bfP_1)
   rn_P0 = -2*diag(BGz1)*E1;%tstr=[rn,rr*[1;P1]] 
   tn_P1 = +2*diag(BGz2)*E2;%tstr=[rn,rr*[1;P1]] 

   %%matrices to map uu to rn & tn coefficients
   M_u2r = 2*diag(BG1)*F1.';
   M_u2t = -2*diag(BG2)*F2.';

   %%matrices to apply the edge conditions;
   M_PM1 = E1.'*diag(-1./Lam1);
   M_PM2 = E2.'*diag(-1./Lam2);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%5th output:
   intrinsic_admittance = ( BG1(1)/BG2(1) );
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else%%SEP_FXN==1
   input1   = {gam1,alp1,H1;
               gam2,alp2,H2};
   input2   = {del0,Dr,sigr,nunu_tilde};
   NMM      = [Nterms,1,0];
   %%
   [MK,forcing,xtra,intrinsic_admittance] =...
      SUB_RTstep_kernel_forcing(input1,input2,NMM,INC_SUB);

   %forcing  = {finc1,ME1;
   %            finc2,ME2};
   finc1 = forcing{1,1};
   ME1   = forcing{1,2};
   finc2 = forcing{2,1};
   ME2   = forcing{2,2};

   %xtra  = {rn_P0,tn_P1;
   %         M_PM1,M_PM2;
   %         M_u2r,M_u2t};
   rn_P0 = xtra{1,1};
   tn_P1 = xtra{1,2};
   M_PM1 = xtra{2,1};
   M_PM2 = xtra{2,2};
   M_u2r = xtra{3,1};
   M_u2t = xtra{3,2};
end

%%SOLVE INTEGRAL EQN:
%uu    = -MK\[1i*F1(:,1),ME1,ME2];
%uu    = -MK\[-1i*finc1,ME1,ME2];
%func  = -1i*[-finc1,ME1,ME2];func(1:10,:)
uu    = -1i*MK\[finc1,ME1,ME2];%uu(1:10,:)
%%
%rr          = 2*diag(BG1)*F1.'*uu;
rr          = M_u2r*uu;
rr(1)       = rr(1)+1;
rr(:,2:3)   = rr(:,2:3)+rn_P0;
   %-2*diag(BGz1)*E1;%tstr=[rn,rr*[1;P1]]
%%
%tt          = -2*diag(BG2)*F2.'*uu;
tt          = M_u2t*uu;
tt(:,4:5)   = tt(:,4:5)+tn_P1;
   %+2*diag(BGz2)*E2;%tstt=[tn tt*[1;P1]]
%%

%%APPLY EDGE CONDITIONS:
if Dr(1)==0%%water on left:
   j_dis        = 2:4;
   rr(:,j_dis)  = [];
   tt(:,j_dis)  = [];
   %%
   M2  = M_PM2(2,:)*tt;
   Q2  = -M2(1)/M2(2);
   %%
   rn  = rr*[1;Q2];
   tn  = tt*[1;Q2];
   %%
   v_unk        = [rn(1);tn(1);Q2];
   uu(:,j_dis)  = [];
   umG          = uu*[1;Q2];
   uuG          = uu;
%  save RT2G v_unk umG uuG
%  tst1=[rn(1),1+2*BG1(1)*F1(:,1).'*umG]
elseif bc==1%%free edges:
   rr(:,[2 4])  = [];
   tt(:,[2 4])  = [];
   %%
   MM     = [M_PM1(2,:)*rr; M_PM2(2,:)*tt];
   MM(1)  = MM(1)+M_PM1(2,1);
   QQ     = -MM(:,2:3)\MM(:,1);
   %%
   rn  = rr*[1;QQ];
   tn  = tt*[1;QQ];
else%%frozen edges:
   rr(:,2:3) = rr(:,2:3)+rr(:,4:5);
   rr(:,4:5) = [];
   %%
   tt(:,2:3) = tt(:,2:3)+tt(:,4:5);
   tt(:,4:5) = [];
   %%
   PM1       = M_PM1*rr;
   PM1(:,1)  = PM1(:,1)+M_PM1(:,1);
   PM2       = M_PM2*tt;
   dPM       = PM2-PM1;
   Q1        = -dPM(:,2:3)\dPM(:,1);
   %%
   rn  = rr*[1;Q1];
   tn  = tt*[1;Q1];
end
  
%%REFLECTION AND TRANSMISSION COEFFICIENTS:
R  = rn(1);
T  = tn(1);
%s  = ( BG1(1)/BG2(1) );%intrinsic admittance
s  = intrinsic_admittance;
y  = {Z_out,rn,tn,s};

if DO_SWAP
   R   = -R'*T/T';
   T   = (1-abs(R)^2)/T';
   s   = 1/s;
   %%
   y(2:3) = y([3 2]);%%swap rn,tn
   y{4}   = s;
end

if do_test

   disp('R&T, |R|&|T|, |R|^2+s*|T|^2:');
   disp([R,T])
   disp([abs(R),abs(T)])
   %disp( abs(gam1(1)-gam2(1))/(gam1(1)+gam2(1)) )
   disp(R*R'+s*T*T')
   %%
   fp1 = gam1.^2+nunu_tilde(1);
   fp2 = gam2.^2+nunu_tilde(2);
   fm1 = gam1.^2-nunu_tilde(1);
   fm2 = gam2.^2-nunu_tilde(2);
   %%
   del1  = lam-sigr(1)*mu;
   Lam1  = Dr(1)*gam1.^4+del1;
   del2  = lam-sigr(2)*mu;
   Lam2  = Dr(2)*gam2.^4+del2;

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
 			phys_vars,hh([2 1]),bc,NN,INC_SUB,EE(:,[2 1]),rho_wtr);
      disp('actual R&T:')
      disp([Rm2 Tm2]);
   end
end

% function y = calc_res(Z2,hr,gamma)
% %% y=calc_res(Z2,hr,gamma)=Res(1/f(K),gamma_n),
% %% where gamma_n is a root of the dispersion relation
% %% f=1/K/tanh(KH)-(Dr*K^4+lam-mr*mu);
% %% Z2={lam,mu,H}.
% lam   = Z2{1};
% mu    = Z2{2};
% H     = Z2{3};
% Dr    = hr^3;
% mr    = hr;
% %%
% Gam   = Dr*gamma.^4+lam-mr*mu;
% Gampr = Gam+4*Dr*gamma.^4;
% denom = H*(Gam.^2.*gamma.^2-1)+Gampr;
% y     = -gamma./denom;

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
