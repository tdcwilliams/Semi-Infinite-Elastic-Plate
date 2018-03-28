%%inputs to SUB_RTstep_Gal_comp.m
period      = 10;%% wave period [s]
theta_inc   = 0;%% wave incident angle [degrees]
H_dim       = 100;%% water depth [m]
phys_vars   = {period,theta_inc,H_dim};

if 1%%semi-inf plate
   hh = [0 1];
   bc = 1;
end

MM = [1,1];
NN = [10,1000];
youngs = 5.45e9;

[out,y] = SUB_RTstep_Gal_comp(...
   phys_vars,hh,bc,NN,MM,youngs);

Aci1  = 0;
Aci2  = 0;
Ainc1 = 1;
Ainc2 = 0;

%% run old program with these incident waves
%% and inputs
Ai_vec2  = [Aci2,Ainc2].';
ACnew = [0,0];
if hh(1)==0
   %%waves from left: [flex-grav amp]
   Ai_vec1 = Ainc1;
   Ai_vec   = [0;Ai_vec1;Ai_vec2];
else
   %%waves from left: [compressive amp,flex-grav amp]
   Ai_vec1 = [Aci1,Ainc1].';
   ACnew(1) = out.Rp(1,:)*Ai_vec1 +...
               out.Tm(1,:)*Ai_vec2;
   Ai_vec   = [Ai_vec1;Ai_vec2];
end

ACnew(2) = out.Tp(1,:)*Ai_vec1 +...
            out.Tm(1,:)*Ai_vec2;
R        = out.Rp(end,:)*Ai_vec1 +...
            out.Tm(end,:)*Ai_vec2;
T        = out.Tp(end,:)*Ai_vec1 +...
            out.Tm(end,:)*Ai_vec2;
RTnew    = [R,T];

%%waves from right: [compressive amp,flex-grav amp]
out0 = SUB_RTstep_Gal_comp_v1(...
         phys_vars,hh,bc,NN,youngs,Ai_vec);
ACold = [out0.Ac_scat,out0.Bc_scat]
ACnew

RTold = [out0.R,out0.T]
RTnew
