%%inputs to SUB_RTstep_Gal_comp.m
period      = 10;%% wave period [s]
theta_inc   = 0;%% wave incident angle [degrees]
H_dim       = 100;%% water depth [m]
phys_vars   = {period,theta_inc,H_dim};

if 1%%semi-inf plate
   hh = [0 1];
   bc = 1;
end

NN       = [10,1000];
youngs   = 5.45e9;

if 1
   [RT_mats,Ac_scat,y]  = SUB_RTstep_Gal_comp(...
                           phys_vars,hh,bc,NN,youngs);

   %%scattered flex-grav waves;
   Rp = RT_mats{1,1};%%1x2 matrix
   Tp = RT_mats{2,1};%%1x2 matrix
   Rm = RT_mats{1,2};%%1x2 matrix
   Tm = RT_mats{2,2};%%1x2 matrix

   %%scattered compressive waves;
   Acs1p = Ac_scat{1,1};
   Acs2p = Ac_scat{2,1};
   Acs1m = Ac_scat{1,2};
   Acs2m = Ac_scat{2,2};
end

Aci1  = 0;
Aci2  = 0;
Ainc1 = 1;
Ainc2 = 0;

%% run old program with these incident waves
%% and inputs
Ai_vec1  = [Aci1,Ainc1].';%%waves from right: [compressive amp,flex-grav amp]
Ai_vec2  = [Aci2,Ainc2].';%%waves from left
Ai_vec   = [Ai_vec1,Ai_vec2];
[R,T]    = SUB_RTstep_Gal_comp_v1(...
               phys_vars,hh,bc,NN,youngs,Ai_vec);
RTold    = [R,T]

%% display results from new program
Acs1     = Acs1p*Ai_vec1 + Acs1m*Ai_vec2
Acs2     = Acs2p*Ai_vec1 + Acs2m*Ai_vec2
R        = Rp*Ai_vec1    + Rm*Ai_vec2
T        = Tp*Ai_vec1    + Tm*Ai_vec2
RTnew    = [R,T]



