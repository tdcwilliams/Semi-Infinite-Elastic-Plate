function [R,T,y] = SUB_RTstep_Gal(...
			phys_vars,hh,bc,NN,INC_SUB,EE,rho_wtr,DO_KC)
%% calc's scattering coefficients for an ice step
%% CALL: [R,T,y]=...
%%   SUB_RTstep_Galerkin(phys_vars,hh,bc,NN,INC_SUB,DO_KC)
%% INPUTS: N is no of imag roots to use;
%%   phys_vars=(vector) T or [T, theta_inc=angle of incidence];
%%     or (cell) {T,theta_inc,H_dim}.
%%   where H_dim is distance from sea floor
%%     to bottom of thickest ice sheet;
%%   hh  = [h1,h2] are ice thicknesses for LHS & RHS
%%   bc  = 0: frozen edge conditions; 1: free edges
%%   NN  = vector[Nterms,Nroots]
%%     -> no of polys to use for P' in Galerkin expansion
%%     & no of roots to use
%%   INC_SUB: use submergence or not
%%   EE = [E1,E2;rho1,rho2;nu1,nu2],
%%    where E*, rho*, nu* are ice Young's modulus, density and Poisson's ratio (1: LHS,2: RHS);
%%    or E (same Young's modulus on both sides) (& default values for YM and rho)
%%   rho_wtr  = water density
%%   DO_KC = 1: do kernel correction; 0: don't
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
else
   if prod(size(EE))==1
      EE = [EE,EE;
            prams(4),prams(4);
            prams(5),prams(5)];
   end
end
if ~exist('rho_wtr')
   rho_wtr  = prams(3);
end
if ~exist('DO_KC')
   DO_KC = 1;
end

do_test  = 0;
if nargin==0
   %% do some tests
   %% - print |R|&|T|, and check energy;
   do_test  = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WANT LARGER ICE THICKNESS ON RIGHT:
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%CALCULATE KERNEL MATRIX AND FORCING VECTORS
input1   = {gam1,alp1,H1;
            gam2,alp2,H2};
input2   = {del0,Dr,sigr,nunu_tilde};
NMM      = [Nterms,1,0];
%%
[MK,forcing,xtra,intrinsic_admittance] =...
   SUB_step_GAL_kernel_forcing(input1,input2,NMM,INC_SUB,DO_KC);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%SOLVE INTEGRAL EQN:
uu = MK\[finc1,ME1,ME2];%uu(1:10,:)
y  = {Z_out,[],[],intrinsic_admittance,uu};
%%
rr          = M_u2r*uu;%rr(1:10,:)
rr(1)       = rr(1)+1;
rr(:,2:3)   = rr(:,2:3)+rn_P0;%rr(1:10,:)
%%
tt          = M_u2t*uu;
tt(:,4:5)   = tt(:,4:5)+tn_P1;

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
s        = intrinsic_admittance;
y(2:3)   = {rn,tn};

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
