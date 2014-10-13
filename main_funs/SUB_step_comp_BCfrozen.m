%% SUB_step_comp_BCwtr.m
%% Author: Timothy Williams
%% Date: 20141010, 16:55:38 CEST
function [ScatMats,y2] =...
   SUB_step_comp_BCfrozen(rr2,tt2,MM,input_bc,SURGE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%book-keeping
M1       = MM(1);
M2       = MM(2);
%%
jd0      = [3 4];%%S(0^+)=S(0^-), \psi(0^+)=\psi(0^-), so remove "+"
Ninc     = M1+M2+2;%%no of incident waves
j_unk    = [Ninc+(1:4)];
j_inc    = 1:Ninc;
j_dis2   = Ninc+jd0;
%%
rr2(:,j_dis2-2)  = rr2(:,j_dis2-2)+rr2(:,j_dis2);
tt2(:,j_dis2-2)  = tt2(:,j_dis2-2)+tt2(:,j_dis2);
rr2(:,j_dis2)  = [];
tt2(:,j_dis2)  = [];
%%
jr_inc1  = 1:M1;        %%rows - f-g waves scattered to  LHS
jc_inc1  = 1:M1;        %%cols - f-g waves incident from LHS
jr_inc2  = 1:M2;        %%rows - f-g waves scattered to  RHS
jc_inc2  = M1+(1:M2);   %%cols - f-g waves incident from RHS
JC_comp  = M1+M2+(1:2); %%cols - compressive wave incident from RHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
kc1   = input_bc{1,1};
Kc1   = input_bc{1,2};
kc2   = input_bc{2,1};
Kc2   = input_bc{2,2};
zc1   = input_bc{3,1};
zc2   = input_bc{3,2};
dzc   = zc1-zc2;
%%
M_PM1B   = input_bc{4,1}.';
M_PM2B   = input_bc{4,2}.';
M_M0w    = input_bc{5,1}.';
M_M1w    = input_bc{5,2}.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Medge = zeros(4,size(rr2,2));

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

if SURGE==0
   if dzc>0
      %% get 2 more eqn's: u(0+)=u(0-)=0
      %% => Ac_inc,Bc_inc are now unknowns
      %% these are now the amplitudes of standing waves:
      %% 2i*Ac_inc*sin(kc1*x),-2i*Bc_inc*sin(kc2*x)

      %%last 2 columns are now multiplied by 0
      jd3            = j_unk(3:4);
      M_edge(:,jd3)  = [];
      rr2(:,jd3)     = [];
      tt2(:,jd3)     = [];

      %%2 fewer inc waves as compressional wave amp's are now unknowns
      j_inc(end-1:end)  = [];
      j_unk             = [JC_comp,j_unk(1:2)];%% 1st 2 unknowns are now Ac_inc,Bc_inc
   else
      %% Get 2 more eqn's: u(0+)=u(0-)=0,
      %% but 1st eqn (u_1 cts) is redundant.
      %% In this case f-g/comp problems are decoupled anyway,
      %% so remove comp problem and solve it analytically.
      %% Still have standing waves:
      %% 2i*Ac_inc*sin(kc1*x),-2i*Bc_inc*sin(kc2*x),
      %% but M0 cts => kc1*Ac_inc=-kc2*Bc_inc

      %%last 2 columns are now multiplied by 0
      jd3            = [JC_comp,j_unk(3:4)];
      M_edge(:,jd3)  = [];
      rr2(:,jd3)     = [];
      tt2(:,jd3)     = [];

      %%2 fewer eqns
      M_edge(1:2,:)     = [];

      %%2 fewer inc waves
      %%2 fewer unknowns
      j_inc(end-1:end)  = [];
      j_unk             = j_unk(1:2);%%[S(0\pm),psi(0\pm)]
   end
end

Q2B   = -Medge(:,j_unk)\Medge(:,j_inc);
rn2   = rr2(:,j_inc);
tn2   = tt2(:,j_inc);
rn2   = rn2+rr2(:,j_unk)*Q2B;
tn2   = tn2+tt2(:,j_unk)*Q2B;
%%
if SURGE==1
   Ac_scat              = Q2B(3,:);%% u(0-) = Ac_inc+Ac_scat
   Ac_scat(JC_comp(1))  = Ac_scat(JC_comp(1))-1;
   Bc_scat              = Q2B(4,:);%% u(0+) = Bc_scat+Bc_inc
   Bc_scat(JC_comp(2))  = Bc_scat(JC_comp(2))-1;

   %%split according to different groups of inc waves;
   Ac_split = {Ac_scat(:,JC_comp(1)),
               Ac_scat(:,jc_inc1),
               Ac_scat(:,JC_comp(2)),
               Ac_scat(:,jc_inc2)};
   rn_split = {rn2(:,JC_comp(1)),
               rn2(:,jc_inc1),
               rn2(:,JC_comp(2)),
               rn2(:,jc_inc2)};
   Bc_split = {Bc_scat(:,JC_comp(1)),
               Bc_scat(:,jc_inc1),
               Bc_scat(:,JC_comp(2)),
               Bc_scat(:,jc_inc2)};
   tn_split = {tn2(:,JC_comp(1)),
               tn2(:,jc_inc1),
               tn2(:,JC_comp(2)),
               tn2(:,jc_inc2)};

   %% unknowns in standard order
   Q2B   = [Q2B(1:2,:); %%S(0-),psi(0-)
            Q2B(3,:);   %%u(0-)
            Q2B(1:2,:); %%S(0+),psi(0+)
            Q2B(4,:)];  %%u(0+)
elseif dzc>0
   %%calc amplitudes of standing waves;
   Ac_stand = 2i*Q2B(1,:);
   Bc_stand = -2i*Q2B(2,:);

   %%split according to different groups of inc waves;
   Ac_split = {Ac_stand(:,jc_inc1),
               Ac_stand(:,jc_inc2)};
   rn_split = {rn2(:,jc_inc1),
               rn2(:,jc_inc2)};
   Bc_split = {Bc_scat(:,jc_inc1),
               Bc_scat(:,jc_inc2)};
   tn_split = {tn2(:,jc_inc1),
               tn2(:,jc_inc2)};

   %% unknowns in standard order
   zz    = 0*Q2B(1,:);
   Q2B   = [Q2B(1:2,:);zz;  %% S(0-),psi(0-),u(0-)=0
            Q2B(1:2,:);zz]; %% S(0+),psi(0+),u(0+)=0
else
   %%calc amplitudes of standing waves;
   Ac_split = 1;
   Bc_split = -kc1/kc2*Ac_split;

   %%split according to different groups of inc waves;
   rn_split = {rn2(:,jc_inc1),
               rn2(:,jc_inc2)};
   tn_split = {tn2(:,jc_inc1),
               tn2(:,jc_inc2)};

   %%unknowns in standard order
   zz    = 0*Q2B(1,:);
   Q2B   = [Q2B;zz;  %% S(0-),psi(0-),u(0-)=0
            Q2B;zz]; %% S(0+),psi(0+),u(0+)=0
end
%%
if SURGE==1
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
else
   Rp = rn2(jr_inc1,jc_inc1);
   Tp = tn2(jr_inc2,jc_inc1);
   %%
   Tm = rn2(jr_inc1,jc_inc2);
   Rm = tn2(jr_inc2,jc_inc2);
end

%% Main output:
ScatMats = {Rp,Tm;
            Tp,Rm};

%% extra outputs - used for testing, plotting displacement etc            
y2 = {Ac_split,rn_split,Bc_split,tn_split,Q2B};
