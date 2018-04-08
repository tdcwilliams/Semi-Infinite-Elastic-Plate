classdef SingleInterfaceCompGal
properties
   lhs_info
   rhs_info
   params
   solution
end
methods
   function obj = SingleInterfaceCompGal(...
      phys_vars,hh,bc,NN,MM,youngs)

      %% constructor
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
         %hh = [1 1];
      end
      if ~exist('bc')
         %bc = 1;%% free edge conditions;
         bc = 0;%% frozen edge conditions;
      end
      if ~exist('NN')
         NN = [10 1000];%% [N_poly, N_roots];
      end
      if ~exist('MM')
         MM = [1 1];%% [M1,M2]
      end
      if ~exist('youngs')
         youngs   = 5.45e9*[1 1];
         %youngs   = 5.45e9*[1 .8];
         %youngs   = [1 1]*NDphyspram(1);
      elseif length(youngs)==1
         youngs   = [1 1]*youngs;
      end

      obj.params.period = phys_vars{1};
      obj.params.theta_inc = phys_vars{2};
      obj.params.H_dim = phys_vars{3};
      obj.params.N_poly = NN(1);
      obj.params.edge_condition = bc;
      obj.params.FAST_KERNEL = 0;
      obj.params.INC_SUB = 1;


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
        obj.params.N_roots = NN(2);
        Ninput = obj.params.N_roots*ones(1,1+obj.params.INC_SUB);
      end
      if obj.params.INC_SUB==1
         [Rts,rts,HH,...
            obj.params.alp_y,obj.params.del0,obj.params.L] =...
               NDmakeZ2_sub(phys_vars,hh,Ninput,youngs);
      else
         [Rts,rts,H,obj.params.alp_y,obj.params.del0,obj.params.L] =...
            NDmakeZ2_rel(phys_vars,hh,Ninput,youngs);
         HH = [H H];
      end

      obj.lhs_info.thickness = hh(1);
      obj.lhs_info.depth_nondim = HH(1);
      obj.lhs_info.youngs = youngs(1);
      obj.lhs_info.num_inc_waves = MM(1);
      obj.lhs_info.roots = Rts{1};
      obj.lhs_info.alp_x = rts{1};

      obj.rhs_info.thickness = hh(2);
      obj.rhs_info.depth_nondim = HH(2);
      obj.rhs_info.youngs = youngs(2);
      obj.rhs_info.num_inc_waves = MM(2);
      obj.rhs_info.roots = Rts{2};
      obj.rhs_info.alp_x = rts{2};

      %%want larger ice thickness on right:
      DO_SWAP = ( hh(1)>hh(2) );
      if DO_SWAP
         tmp = obj.lhs_info;
         obj.lhs_info = obj.rhs_info;
         obj.rhs_info = tmp;
         clear tmp;
      end
      obj.params.h_ref = obj.rhs_info.thickness;
      obj.params.E_ref = obj.rhs_info.youngs;
      obj.params.N_roots = length(obj.rhs_info.roots)-3;%number of imaginary roots

      %% common to both sides
      obj.params.nu = NDphyspram(5);
      obj.params.rho_wtr = NDphyspram(3);
      obj.params.rho_ice = NDphyspram(4);
      obj.params.omega = 2*pi/obj.params.period;
      obj.params.lam = obj.params.del0(1);
      obj.params.mu = -obj.params.del0(2);
      obj.params.nu1 = (1-obj.params.nu)*obj.params.alp_y^2;

      %% get order of singularity
      obj.params.mvec = (0:obj.params.N_poly)';
      obj.params.alpC = .5-1/3*obj.params.INC_SUB*(hh(1)~=hh(2));%%=1/6 if submergence included;
			       %%=1/2 if no submergence.

      %% get extend lhs and rhs info
      obj.lhs_info = obj.extend_info(obj.lhs_info);
      obj.rhs_info = obj.extend_info(obj.rhs_info);
   end%% end constructor

   function info = extend_info(obj,info,hmax)

      %%Compressional stuff for u problem:
      info.mu_lame = info.youngs/2/(1+obj.params.nu);
      info.Kc_dim  = info.youngs*info.thickness/(1-obj.params.nu^2);
         %%compressional rigidity ~ Pa*m ~ rho_wtr*om^2*L^2*L
      info.Kc = info.Kc_dim/(obj.params.rho_wtr*obj.params.omega^2*obj.params.L^2)/obj.params.L;
      info.m_dim   = obj.params.rho_ice*info.thickness;
      info.sig_dim = info.m_dim/obj.params.rho_wtr;
      info.zc_dim  = -info.sig_dim+info.thickness/2;

      %% non-dimensionalise
      info.sig = info.sig_dim/obj.params.L;
      info.zc = info.zc_dim/obj.params.L;

      %% wave numbers for compressional problem
      if info.Kc==0
         info.kc_dim = 0;
      else
         info.kc_dim = sqrt(obj.params.omega^2*info.m_dim/info.Kc_dim);
      end
      info.kc = info.kc_dim*obj.params.L;

      %% rotational inertia
      info.J_dim = obj.params.rho_ice*info.thickness^3/12;

      %%make real parts of gam be >0 (for M0_wtr,M1_wtr)
      jneg = find(real(info.roots)<0);
      info.roots(jneg) = -info.roots(jneg);

      info.hr = info.thickness/obj.params.h_ref;
      info.Dr = info.youngs/obj.params.E_ref*info.hr^3;

      %% coefficients in expansion of Green's function's derivatives
      info.BGzz = obj.calc_res(info).*info.roots./info.alp_x;%%BG_0^(0)
      info.Lam  = info.Dr*info.roots.^4+obj.params.lam-info.hr*obj.params.mu;
      info.BGz  = info.Lam.*info.BGzz;
      info.BG   = info.Lam.*info.BGz;%%BG_0^(2)

      %% inner products of polynomials with eigenfunctions
      kap = -1i*info.roots*info.depth_nondim;
      [KAP,MU] = meshgrid(kap,2*obj.params.mvec+obj.params.alpC);
      besJ = besselj(MU,KAP);%% length(mvec) x length(kap1)
      c_left = gamma(obj.params.alpC)*(obj.params.alpC+2*obj.params.mvec).*(-1).^obj.params.mvec;
      c_right = (2./kap).^obj.params.alpC./cos(kap);
      info.F = diag(c_left)*besJ*diag(c_right);%%same as (37)

      %% for edge conditions
      info.fm = info.roots.^2-obj.params.nu1;
      info.fp = info.roots.^2+obj.params.nu1;
      info.E = [-1+0*info.roots,info.Dr*info.fm];
         %% S=-(L^2/D1)*D*w_xxx=D/D1*\pa_{xND}^3(wND), psi=-w_x
      info.ME = info.F*diag(1i*info.BGz)*info.E;
      info.M_PM = diag([-1,1])*info.E.'*diag(1./info.Lam);
         %% M=-L*Dr*Lm(w)=-Dr*-fm*(wND)=Dr*fm=+L*Dr*fm

      %% for energy test
      info.Ew_fac = -.5/info.BG(1);
      info.Ec_fac = info.kc*info.Kc;
   end

   function obj = assemble(obj)
      %% further common parameters
      %dzc_dim  = obj.lhs_info.zc_dim+...
                        -obj.rhs_info.zc_dim;%pause
      %dzc = dzc_dim/obj.params.L;
      obj.params.H = obj.rhs_info.depth_nondim+obj.rhs_info.sig;

      obj.params.M1_rel = obj.params.rho_wtr*obj.params.omega^2*obj.params.L^4;
         %%scale bending moment by this, since
         %%\int\sig_ij.(z-z_c).dz ~ Pa*m^2 ~ rho_wtr*om^2*L^2*L^2 = D1/L
      obj.params.M0_rel = obj.params.rho_wtr*obj.params.omega^2*obj.params.L^3;
         %%scale zero-th moment by this, since
         %% \int\sig_ij.dz ~ Pa*m ~ rho_wtr*om^2*L^2*L = D1/L^2

      %%MAIN KERNEL MATRIX:
      if 1%~FAST_KERNEL%%basic way:-other ways don't work
         obj.solution.MK =...
            obj.lhs_info.F*diag(1i*obj.lhs_info.BG)*obj.lhs_info.F.'+...
               + obj.rhs_info.F*diag(1i*obj.rhs_info.BG)*obj.rhs_info.F.';
               %%looks same as (38c)
      end

      %%Moments of pressure M0w and M1w:
      %%\int_{-\sig2}^{-sig1}\varf
      gam_H_lhs = obj.lhs_info.roots*obj.lhs_info.depth_nondim;
      gam_H_rhs = obj.lhs_info.roots*obj.rhs_info.depth_nondim;
      denom    = 1+exp(-2*gam_H_lhs);
      num0     = 1-exp(-2*gam_H_lhs);
      exfac    = exp(gam_H_rhs-gam_H_lhs);
      M0_wtr = 1./obj.lhs_info.roots.^2./obj.lhs_info.Lam+...
                  -1./obj.lhs_info.roots.*exfac.*num0./denom;
      %%
      M1a = (obj.lhs_info.depth_nondim./obj.lhs_info.Lam-1)./obj.lhs_info.roots.^2;
      num1 = (1-gam_H_rhs)+(1+gam_H_rhs).*exp(-2*gam_H_rhs);
      M1b = num1.*exfac./obj.lhs_info.roots.^2./denom;
      M1_wtr = M1a+M1b-(obj.params.H+obj.rhs_info.zc)*M0_wtr;

      %%contribution from side integrals to rr2:
      obj.solution.rr_wtr_psi2 = -2i*obj.lhs_info.BG.*M1_wtr;
         %%Q2 (psi(0+)) column
      obj.solution.rr_wtr_U2 = -2i*obj.lhs_info.BG.*M0_wtr;
         %%U2 (u(0+)) column

      %% book-keeping
      M1 = obj.lhs_info.num_inc_waves;
      M2 = obj.rhs_info.num_inc_waves;
      Minc_tot = M1+M2+2;
      %%
      obj.solution.j_inc2  = (1:M2); %%col's for inc water waves from lhs
      obj.solution.jc_inc2 = M2+1; %%col for inc compressive wave from lhs
      obj.solution.j_inc1  = M2+1+(1:M1); %%col's for inc water waves from rhs
      obj.solution.jc_inc1 = M2+M1+2; %%col for inc compressive wave from rhs
      %%
      obj.solution.j_Q2 = Minc_tot+(1:2); %%col's for [S(0-),psi(0-)] TODO check which is which
      obj.solution.j_U2 = Minc_tot+3;     %%col for u(0-)
      obj.solution.j_Q1  = Minc_tot+(4:5); %%col's for [S(0+),psi(0+)]
      obj.solution.j_U1  = Minc_tot+6;     %%col for u(0+)
      obj.solution.j_side_ints = [obj.solution.j_Q2(2),...
         obj.solution.j_U2];%%col's corresponding to psi(0+), u(0+)

      %%matrix M0w=M_M0w*rr2
      obj.solution.M_M0w = -M0_wtr.';%%\int\sig_11 dz=\int(-P)dz=-\int\phi dz
      obj.solution.M_M1w = -M1_wtr.';%%\int\sig_11*(z-zc1)dz=\int(-P)*(z-zc1)dz=-\int\phi*(z-zc1)dz

      %%extra forcing from side integral terms in rr2;
      %% col's of forcing_u correspond to
      %% [Ainc2,Ac_inc,Binc2,Bc_inc,...
      %    psi(0-),S(0-),psi(0+),S(0+),u(0-),u(0+)];
      ZM = 0*obj.lhs_info.F(:,1);
      obj.solution.forcing_u = [obj.lhs_info.F(:,1:obj.lhs_info.num_inc_waves),...
                   ZM,...
                   -obj.rhs_info.F(:,1:obj.rhs_info.num_inc_waves),...
                   ZM,...
                   obj.lhs_info.ME,...
                   obj.rhs_info.ME,...
                   ZM,ZM];
      obj.solution.forcing_u(:,obj.solution.j_side_ints) =...
         obj.solution.forcing_u(:,obj.solution.j_side_ints)+...
            +obj.lhs_info.F*[obj.solution.rr_wtr_psi2,...
                              obj.solution.rr_wtr_U2];
   end%% end assemble

   function obj = solve(obj)

      %% assemble matrices and forcing vectors
      obj = obj.assemble();

      %% solve the integral eqn
      obj.solution.u = obj.solution.MK\obj.solution.forcing_u;

      % ================================================
      %%scattered waves for LHS
      obj.solution.lhs.fluid_coeffs =...
         -2i*diag(obj.lhs_info.BG)*...
         obj.lhs_info.F.'*obj.solution.u;%rr2

      %%contrib from inc water waves
      j_inc1 = obj.solution.j_inc1;
      obj.solution.lhs.fluid_coeffs(j_inc1,j_inc1) =...
         + obj.solution.lhs.fluid_coeffs(j_inc1,j_inc1)+...
         + eye(obj.lhs_info.num_inc_waves);

      %%contrib from Q1
      %%**CHANGE FOR MINDLIN**
      j_Q1 = obj.solution.j_Q1;
      obj.solution.lhs.fluid_coeffs(:,j_Q1) =...
         + obj.solution.lhs.fluid_coeffs(:,j_Q1)+...
         + 2i*diag(obj.lhs_info.BGz)*obj.lhs_info.E;

      %%from side integrals;
      obj.solution.lhs.fluid_coeffs(:,obj.solution.j_side_ints) =...
         + obj.solution.lhs.fluid_coeffs(:,obj.solution.j_side_ints)+...
         + [obj.solution.rr_wtr_psi2,obj.solution.rr_wtr_U2];
      % ================================================

      % ================================================
      %%scattered waves for RHS
      obj.solution.rhs.fluid_coeffs =...
         -2i*diag(obj.rhs_info.BG)*...
         obj.rhs_info.F.'*obj.solution.u;%rr2

      %%contrib from inc water waves
      j_inc2 = obj.solution.j_inc2;
      obj.solution.rhs.fluid_coeffs(j_inc2,j_inc2) =...
         + obj.solution.rhs.fluid_coeffs(j_inc2,j_inc2)+...
         + eye(obj.rhs_info.num_inc_waves);

      %%contrib from Q2
      %%**CHANGE FOR MINDLIN**
      j_Q2 = obj.solution.j_Q2;
      obj.solution.rhs.fluid_coeffs(:,j_Q2) =...
         + obj.solution.rhs.fluid_coeffs(:,j_Q2)+...
         - 2i*diag(obj.rhs_info.BGz)*obj.rhs_info.E;
      % ================================================

      if obj.lhs_info.Dr==0
         obj = obj.edge_cons_ice_water();
      elseif obj.params.edge_condition == 0
         error(['Unimplemented edge condition option:',...
            num2str(obj.params.edge_condition)]);
         obj = obj.edge_cons_frozen();
      elseif obj.params.edge_condition == 1
         error(['Unimplemented edge condition option:',...
            num2str(obj.params.edge_condition)]);
         obj = obj.edge_cons_free();
      else
         error(['Unknown edge condition option:',...
            num2str(obj.params.edge_condition)]);
      end

   end%% end solve

   function obj = edge_cons_ice_water(obj)
      j_inc_all = [obj.solution.j_inc1,...
                  obj.solution.j_inc2,...
                  obj.solution.jc_inc2];
         %% no more inc compressive wave from left
      j_unknown = [obj.solution.j_Q2(2),...
                  obj.solution.j_U2];
         %%keep psi(0+), and u(0+)
      M1 = obj.lhs_info.num_inc_waves;
      M2 = obj.rhs_info.num_inc_waves;
      j_inc1 = obj.solution.j_inc1;
      j_inc2 = obj.solution.j_inc2;
      jc_inc2 = obj.solution.jc_inc2;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%bending moment eqn
      MB = obj.rhs_info.M_PM(2,:)*obj.solution.rhs.fluid_coeffs;%% 1 x length(kap2)

      %%inc water waves from RHS
      MB(j_inc2) = MB(j_inc2)+obj.rhs_info.M_PM(2,1:M2);
         %%'+' since M has even derivatives

      %%M_1^w:
      M1w = obj.solution.M_M1w*obj.solution.lhs.fluid_coeffs;

      %%inc water waves from LHS
      M1w(j_inc1) = M1w(j_inc1)+obj.solution.M_M1w(1:M1);

      Medge = MB-M1w;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%compression at edge:
      %%U2=(Bc_inc+Bc_scat)=(M_0^w/1i/kc2/Kc2)+2*Bc_inc
      M0w = obj.solution.M_M0w*obj.solution.lhs.fluid_coeffs;

      %%inc water waves from LHS
      M0w(j_inc1) = M0w(j_inc1)+diag(obj.solution.M_M0w(1:M1));
      %%
      Bc_inc     = 0;
      Medge(2,:) = M0w/(1i*obj.rhs_info.kc*obj.rhs_info.Kc);

      %%inc compressive wave from RHS
      Medge(2,jc_inc2) = Medge(2,jc_inc2)+2;

      %%adjust u(0+) column
      Medge(2,obj.solution.j_U2) = Medge(2,obj.solution.j_U2)-1;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      obj.solution.QU2 = -Medge(:,j_unknown)\Medge(:,j_inc_all);%%2xNinc matrix

      %% col's of obj.solution.[lhs,rhs].fluid_coeffs tn2 corresp to
      %% [Aw_inc2,Ac_inc2,Aw_inc1]
      %% where incident flex-grav waves are:
      %%  Aw_inc1  = [Ainc1(1),...,Ainc1(M1)]
      %%  Aw_inc2  = [Ainc2(1),...,Ainc2(M2)]
      %% and incident comp wave from rhs is
      %%  Ac_inc2
      %% NB no compressive wave from left (water)
      obj.solution.lhs.fluid_coeffs =...
         obj.solution.lhs.fluid_coeffs(:,j_inc_all)+...
            obj.solution.lhs.fluid_coeffs(:,j_unknown)*...
               obj.solution.QU2;
      obj.solution.rhs.fluid_coeffs =...
         obj.solution.rhs.fluid_coeffs(:,j_inc_all)+...
            obj.solution.rhs.fluid_coeffs(:,j_unknown)*...
               obj.solution.QU2;

      %% compressive wave
      obj.solution.rhs.comp_coeff = obj.solution.QU2(2,:);
         %% - incident compressive wave
      obj.solution.rhs.comp_coeff(jc_inc2) = obj.solution.rhs.comp_coeff(jc_inc2)-1;
      %%
      obj.solution.lhs.comp_coeff = 0*obj.solution.QU2(2,:);
      %%
      obj.solution.u =...
         obj.solution.u(:,j_inc_all)+...
            obj.solution.u(:,j_unknown)*obj.solution.QU2;
   end%%edge_cons_ice_water

   function obj = edge_cons_free(obj)
      j_inc_all = [obj.solution.j_inc1,...
                  obj.solution.j_inc1,...
                  obj.solution.j_inc2,...
                  obj.solution.jc_inc2];
         %% no more inc compressive wave from left
      j_unknown = [obj.solution.j_Q2(2),...
                   obj.solution.j_U2
                   obj.solution.j_Q1(2),...
                   obj.solution.j_U1];
         %%keep psi(0\pm), and u(0\pm)
      M1 = obj.lhs_info.num_inc_waves;
      M2 = obj.rhs_info.num_inc_waves;
      j_inc1 = obj.solution.j_inc1;
      j_inc2 = obj.solution.j_inc2;
      jc_inc1 = obj.solution.jc_inc1;
      jc_inc2 = obj.solution.jc_inc2;

      %% ==============================================================
      %%Bending moment eqn (LHS)
      Medge = obj.lhs_info.M_PM(2,:)*obj.solution.lhs.fluid_coeffs;

      %%inc water waves from LHS
      Medge(1,j_inc1) = Medge(1,j_inc1)+obj.lhs_info.M_PM(2,1:M1);
      %% ==============================================================

      %% ==============================================================
      %%Bending moment eqn (RHS)
      Medge(2,:) = obj.rhs_info.M_PM(2,:)*obj.solution.rhs.fluid_coeffs;

      %%inc water waves from RHS
      Medge(2,j_inc2) = Medge(2,j_inc2)+obj.rhs_info.M_PM(2,1:M2);
         %%'+' since M has even derivatives

      %%M_1^w:
      M1w = obj.solution.M_M1w*obj.solution.lhs.fluid_coeffs;

      %%inc water waves from LHS
      M1w(j_inc1) = M1w(j_inc1)+obj.solution.M_M1w(1:M1);
      Medge(2,:) = Medge(2,:) - M1w;
      %% ==============================================================
   
      %% ==============================================================
      %%compression at edge (LHS):
      %% U1=Ac_inc+Ac_scat=2*Ac_inc
      M_edge(3,:) = 0;

      %% compressive wave from left
      M_edge(3,jc_inc1) = 2;
      M_edge(3,j_U1) = -1;
      %% ==============================================================
   
      %% ==============================================================
      %%compression at edge (RHS):
      %% U1=Ac_inc+Ac_scat=2*Ac_inc
      M_edge(4,:) = 0;

      %% compressive wave from right
      fac = 1i*obj.rhs_info.kc*obj.rhs_info.Kc;
      M_edge(4,jc_inc2) = 2*fac;
      M_edge(4,j_U2) = -1*fac;

      %%inc water waves from LHS
      M0w = obj.solution.M_M0w*obj.solution.lhs.fluid_coeffs;
      M0w(j_inc1) = M0w(j_inc1)+ obj.solution.M_M0w(1:M1);
      M_edge(4,:) = M_edge(4,:) - M0w;
      %% ==============================================================

   end %% edge_cons_free

   function y = calc_res(obj,info)
      %% y=obj.calc_res(Z2,Dr,hr,gamma)=Res(1/f(K),gamma_n),
      %% where gamma_n is a root of the dispersion relation
      %% f=1/K/tanh(KH)-(Dr*K^4+lam-mr*mu);
      %% Z2={lam,mu,H}.
      lam = obj.params.lam;
      mu  = obj.params.mu;
      mr  = info.hr;
      Dr  = info.Dr;
      H = info.depth_nondim;
      gamma = info.roots;
      %%
      Gam   = Dr*gamma.^4+lam-mr*mu;
      Gampr = Gam+4*Dr*gamma.^4;
      denom = H*(Gam.^2.*gamma.^2-1)+Gampr;
      y     = -gamma./denom;
   end%% calc_res

   function [S,Efacs]=get_scattering_matrix(obj)
      j_inc2 = obj.solution.j_inc2;
      jc_inc2 = obj.solution.jc_inc2;
      j_inc1 = obj.solution.j_inc1;

      Efacs = [obj.rhs_info.Ew_fac,...
               obj.rhs_info.Ec_fac,...
               obj.lhs_info.Ew_fac];

      %% water wave from right
      S(:,1) = [obj.solution.rhs.fluid_coeffs(1,j_inc2(1));... %reflected water wave
                obj.solution.rhs.comp_coeff(1,j_inc2(1));...   %generated compressive wave
                obj.solution.lhs.fluid_coeffs(1,j_inc2(1))];   %transmitted water wave

      %% compressive wave from right
      S(:,2) = [obj.solution.rhs.fluid_coeffs(1,jc_inc2(1));... %generated water wave (rhs)
                obj.solution.rhs.comp_coeff(1,jc_inc2(1));...   %reflected compressive wave
                obj.solution.lhs.fluid_coeffs(1,jc_inc2(1))];   %generated water wave (lhs)

      %% water wave from left
      S(:,3) = [obj.solution.rhs.fluid_coeffs(1,j_inc1(1));... %transmitted water wave
                obj.solution.rhs.comp_coeff(1,j_inc1(1));...   %generated compressive wave
                obj.solution.lhs.fluid_coeffs(1,j_inc1(1))];   %reflected water wave

      if obj.lhs_info.hr>0
         jc_inc1 = obj.solution.jc_inc1;

         % generated compressive wave (lhs)
         S(4,:) = [obj.solution.rhs.comp_coeff(1,j_inc2(1));...  %by ww from right
                   obj.solution.rhs.comp_coeff(1,jc_inc2(1));... %by cw from right
                   obj.solution.rhs.comp_coeff(1,j_inc1(1))];    %by ww from left

         % compressive wave from left
         S(:,4) = [obj.solution.rhs.fluid_coeffs(1,jc_inc1(1));... %generated water wave (rhs)
                   obj.solution.rhs.comp_coeff(1,jc_inc1(1));...   %transmitted compressive wave
                   obj.solution.lhs.fluid_coeffs(1,jc_inc1(1));... %generated water wave (lhs)
                   obj.solution.lhs.comp_coeff(1,jc_inc1(1))];     %reflected compressive wave

         Efacs(4) = obj.lhs_info.Ec_fac;
      end
   end%%get_scattering_matrix

   function Etest = test_energy(obj)
      [S,Efacs] = obj.get_scattering_matrix()
      Etest = diag(Efacs)*S*S'*diag(1./Efacs);

      E_in = Efacs(3)
      E_out = Efacs(3)*abs(S(3,3))^2+...
               Efacs(2)*abs(S(2,3))^2+...
               Efacs(1)*abs(S(1,3))^2
   end%% test_energy

end%%end methods
end%%class
