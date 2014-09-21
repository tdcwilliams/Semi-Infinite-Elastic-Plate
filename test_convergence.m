%% test_convergence.m
%% Author: Timothy Williams
%% Date:   20140921, 03:16:04 CEST


period      = 10;%% wave period [s]
theta_inc   = 0;%% wave incident angle [degrees]
H_dim       = 100;%% water depth [m]
phys_vars   = {period,theta_inc,H_dim};
%%
hh = [0 1];
bc = 1;%% free edge conditions;
%bc = 0;%% frozen edge conditions;

INC_SUB  = 1;%% include submergence or not;
prams    = NDphyspram(0);%[E,g,rho_wtr,rho_ice,nu];
EE       = [prams(1),prams(1);
            prams(4),prams(4);
            prams(5),prams(5)];
rho_wtr  = prams(3);

%%% test without kernel correction:
%DO_KC    = 0
%Npolys   = 50;
%Nvec0    = [1000:1000:5000,7e3]';
%for j=1:length(Nvec0)
%   NN                   = [Npolys,Nvec0(j)]
%   [R0(j,1),T0(j,1),y]  = SUB_RTstep_Galerkin(...
%                            phys_vars,hh,bc,NN,INC_SUB,EE,rho_wtr,DO_KC);
%end

%% test with kernel correction:
DO_KC    = 1
Npolys   = 50;
Nvec1    = [250,500:500:2000,5000,7e3]';
for j=1:length(Nvec1)
   NN                   = [Npolys,Nvec1(j)]
   [R1(j,1),T1(j,1),y]  = SUB_RTstep_Galerkin(...
                            phys_vars,hh,bc,NN,INC_SUB,EE,rho_wtr,DO_KC);
end

[R0,T0]
[R1,T1]
