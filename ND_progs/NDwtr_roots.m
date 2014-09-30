function [k_cell,outputs]=NDwtr_roots(inputs,N,type)

%% NON-DIMENSIONALISE:
if type==0
  %% inputs=[gamma,H], dispersion relation is
  %% F(k)=gamma*k*tanh(k*H)-1=0;
  %% put k=K/H
  %%  => (gamma/H)*K*tanh(K)-1=0;
  gam=inputs(:,1);
  H=inputs(:,2);
  gam_on_H=gam./H;
else
  %% inputs=[period,H];
  period=inputs(:,1);
  H=inputs(:,2);
  %%
  om=2*pi./period;
  pram=NDphyspram(0);
  g=pram(2); %m/s^2
  gam_on_H=g./(H.*om.^2)
end

%% CALCULATE OUTPUTS:
M=length(gam_on_H);
k_cell=cell(M,1);
for j=1:M
  k_cell{j}=RTS_wtr_roots(gam_on_H(j),1,N)/H(j);
end
outputs=[gam_on_H,H];