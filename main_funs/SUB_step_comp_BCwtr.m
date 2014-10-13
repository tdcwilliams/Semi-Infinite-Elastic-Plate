%% SUB_step_comp_BCwtr.m
%% Author: Timothy Williams
%% Date: 20141010, 16:55:38 CEST
function [ScatMats,y2] =...
   SUB_step_comp_BCwtr(rr2,tt2,MM,input_bc,SURGE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%book-keeping
M1       = MM(1);
M2       = MM(2);
%%
Ninc     = M1+M2+1;%%no of incident waves
j_inc    = 1:Ninc;
j_unk    = [Ninc+(1:2)];
jr_inc1  = 1:M1;        %%rows - f-g waves scattered to  LHS
jc_inc1  = 1:M1;        %%cols - f-g waves incident from LHS
jr_inc2  = 1:M2;        %%rows - f-g waves scattered to  RHS
jc_inc2  = M1+(1:M2);   %%cols - f-g waves incident from RHS
JC_comp  = M1+M2+1;     %%col  - compressive wave incident from RHS

jd0      = [1 2 3 5];%%no more [S1 (1),Q1 (2),U1 (5)], & S2=0 (3)
                     %%keep Q2 (4) and U2 (6)
jd1      = Ninc;%%no longer have incident compressive wave from left
j_dis2   = [jd1,M1+M2+2+jd0];
%%
rr2(:,j_dis2)  = [];
tt2(:,j_dis2)  = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%kc1   = input_bc{1,1};
%Kc1   = input_bc{1,2};
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

%%split according to different groups of inc waves;
rn_split = {rn2(:,JC_comp),
            rn2(:,jc_inc1),
            rn2(:,jc_inc2)};
tn_split = {tn2(:,JC_comp),
            tn2(:,jc_inc1),
            tn2(:,jc_inc2)};
Bc_split = {Bc_scat(:,JC_comp),
            Bc_scat(:,jc_inc1),
            Bc_scat(:,jc_inc2)};

Rp = rn2(jr_inc1,jc_inc1);
Tp = [Bc_scat(jc_inc1);
      tn2(jr_inc2,jc_inc1)];
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
%else
%   {'R',Rp*Ainc2+Tm(1,:)*[Bc_inc;Binc2];
%    'T',Tp(2,:)*Ainc2+Rm(2,:)*[Bc_inc;Binc2];
%    'Bc_scat',Tp(1,:)*Ainc2+Rm(1,:)*[Bc_inc;Binc2]}
end

%% Main output:
ScatMats = {Rp,Tm;
            Tp,Rm};

%% extra outputs - used for testing, plotting displacement etc            
y2 = {rn_split,Bc_split,tn_split,Q2B};
