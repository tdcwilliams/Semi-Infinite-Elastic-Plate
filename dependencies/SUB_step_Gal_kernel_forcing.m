%% SUB_step_Gal_kernel_forcing.m
%% Author: Timothy Williams
%% Date: 20140915, 10:42:15 CEST
function [MK,forcing,xtra,intrinsic_admittance]  =...
   SUB_step_Gal_kernel_forcing(input1,input2,NMM,INC_SUB,CORRECT_KERNEL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Sort out inputs;

%CORRECT_KERNEL = 0;%%use polylogarithm/zeta function correction to kernel matrix

%input1   = {gam1,alp1,H1;
%            gam2,alp2,H2};
gam1     = input1{1,1};
alp1     = input1{1,2};
H1       = input1{1,3};
gam2     = input1{2,1};
alp2     = input1{2,2};
H2       = input1{2,3};
Nroots   = length(gam2)-3;

%input2   = {del0,Dr,sigr,nunu_tilde};
lam         = input2{1}(1);
mu          = -input2{1}(2);
Dr          = input2{2};
sigr        = input2{3};
nunu_tilde  = input2{4};

%NMM  = [Npolys,M1,M2];
Npolys   = NMM(1);
M1       = NMM(2);
M2       = NMM(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Coefficients of Green's fxn expansions
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

%%intrinsic admittance (for energy conservation check);
intrinsic_admittance = ( BG1(1)/BG2(1) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALC Fm0 & Fmr:
mvec  = (0:Npolys)';
alpC  = .5-1/3*INC_SUB*(sigr(1)~=1);%%=1/6 if submergence included;
			            %%=1/2 if no submergence.

kap1     = -1i*gam1*H2;
kap2     = -1i*gam2*H2;
c_left   = gamma(alpC)*(alpC+2*mvec).*(-1).^mvec;

if 0%%older version of besselj
   besJ1 = besselj(2*mvec'+alpC,kap1).';
   besJ2 = besselj(2*mvec'+alpC,kap2).';
else
   [NU1,Z1] = meshgrid(2*mvec'+alpC, kap1);
   besJ1    = besselj(NU1,Z1).';
   [NU1,Z1] = meshgrid(2*mvec'+alpC, kap2);
   besJ2    = besselj(NU1,Z1).';
end

c_rt1 = (2./kap1).^alpC./cosh(gam1*H1);
F1    = diag(c_left)*besJ1*diag(c_rt1);
c_rt2 = (2./kap2).^alpC./cos(kap2);
F2    = diag(c_left)*besJ2*diag(c_rt2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MAIN KERNEL MATRIX:
MK = F2*diag(BG2)*F2.'+F1*diag(BG1)*F1.';
if CORRECT_KERNEL
   %%try to accelerate convergence of kernel matrix;
   jtest = Nroots-10:Nroots;
   mtest = 2;
   nnvec = (1:Nroots)';

   %%approximate roots;
   gam1ap   = 1i*nnvec*pi/H1;
   gam2ap   = 1i*nnvec*pi/H2;
   %tst_gam1 = [gam1ap(jtest),gam1(jtest+end-Nroots)]
   %tst_gam2 = [gam2ap(jtest),gam2(jtest+end-Nroots)]

   %%approximate BG1,BG2
   BGap   = 1i/pi./nnvec;
   %tst_BG = [BGap(jtest),BG1(jtest+end-Nroots),BG2(jtest+end-Nroots)]

   %%approximate c_left,c_rt1,c_rt2
   Hr    = H1/H2;
   clap  = gamma(alpC)*(alpC+2*mvec).*(-1).^mvec;  %%~c_left: m
   crap1 = (2./(nnvec*pi*Hr)).^alpC.*(-1).^nnvec;  %%~c_rt1:  n
   crap2 = (2./(pi*nnvec)).^alpC.*(-1).^nnvec;     %%~c_rt2:  n

   %%approx bessel functions:
   bet1        = pi/2*(alpC+.5);
   besJ1ap_m   = (-1).^mvec*sqrt(2*Hr/pi^2);
   besJ1ap_n   = nnvec.^(-.5).*cos(nnvec*pi/Hr-bet1);
   besJ1ap     = besJ1ap_m*besJ1ap_n.';
   %tst_besJ1ap = [besJ1ap(mtest,jtest).',besJ1(mtest,jtest+end-Nroots).']
   %%
   besJ2ap_m   = (-1).^mvec*sqrt(2/pi^2);
   besJ2ap_n   = nnvec.^(-.5).*cos(nnvec*pi-bet1);
   besJ2ap     = besJ2ap_m*besJ2ap_n.';
   %tst_besJ2ap = [besJ2ap(mtest,jtest).',besJ2(mtest,jtest+end-Nroots).']

   %%approx F1,F2
   f1m_coeffs  = gamma(alpC)*(alpC+2*mvec)*sqrt(2*Hr/pi^2)*(2*Hr/pi)^alpC;
   f1n_coeffs  = (-1).^nnvec./nnvec.^(alpC+.5).*cos(nnvec*pi/Hr-bet1);
   F1ap        = f1m_coeffs*f1n_coeffs.';
   %tst_F1ap    = [F1ap(mtest,jtest).',F1(mtest,jtest+end-Nroots).']
   f2m_coeffs  = gamma(alpC)*(alpC+2*mvec)*sqrt(2/pi^2)*(2/pi)^alpC;
   f2n_coeffs  = (-1).^nnvec./nnvec.^(alpC+.5).*cos(nnvec*pi-bet1);
   F2ap        = f2m_coeffs*f2n_coeffs.';
   %tst_F2ap    = [F2ap(mtest,jtest).',F2(mtest,jtest+end-Nroots).']

   MKap_N   = F2ap*diag(BGap)*F2ap.'+F1ap*diag(BGap)*F1ap.';
   %%
   ex0   = exp(-2i*bet1);
   ex1   = exp(2i*pi/Hr);
   s_MK  = 2*(alpC+1);
   sfac1 = zeta(s_MK)+real( ex0*SF_polylog(ex1,s_MK) );
   sfac2 = zeta(s_MK)*( 1+real(ex0) );
   MKap  = f2m_coeffs*f2m_coeffs.'*sfac2+...
           f1m_coeffs*f1m_coeffs.'*sfac1; 
   MK    = MK - MKap_N +MKap;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%FORCING TERMS:
fm1   = gam1.^2-nunu_tilde(1);
E1    = [1+0*gam1,-Dr(1)*fm1];
fm2   = gam2.^2-nunu_tilde(2);
E2    = [1+0*gam2,-Dr(2)*fm2];
%%
jinc1 = 1:M1;
jinc2 = 1:M2;
jinc  = 1:(M1+M2);
finc1 = -F1(:,jinc1);
finc2 = F2(:,jinc2);
ME2   = F2*diag(BGz2)*E2;
ME1   = F1*diag(BGz1)*E1;
%%
forcing  = {finc1,ME1;
            finc2,ME2};

%%scattering caused by "forced oscillations" (\bfP_0 & \bfP_1)
rn_P0 = -2*diag(BGz1)*E1;%tstr=[rn,rr*[1;P1]] 
tn_P1 = +2*diag(BGz2)*E2;%tstr=[rn,rr*[1;P1]] 

%%matrices to apply edge conditions:
M_PM1 = E1.'*diag(-1./Lam1);
M_PM2 = E2.'*diag(-1./Lam2);

%%matrices to map uu to rn & tn coefficients
M_u2r = 2*diag(BG1)*F1.';
M_u2t = -2*diag(BG2)*F2.';

xtra  = {rn_P0,tn_P1;
         M_PM1,M_PM2;
         M_u2r,M_u2t};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=calc_res(Z2,gamma)
%% Calculate residue
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
