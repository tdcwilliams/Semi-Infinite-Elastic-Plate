function [k_cell,outputs]=NDice_roots(inputs,N,type)

jw=find(inputs(:,1)==0);
jnz=find(inputs(:,1)~=0);
%%
outputs=0*inputs;
k_cell=cell( size(outputs,1),1 );
if ~isempty(jw)
  [k_cell(jw),outputs(jw,2:3)]=NDwtr_roots(inputs(jw,2:3),N,type);
  inputs(jw,:)=[];
end
if isempty(jnz)
  return;
end

%% NON-DIMENSIONALISE:
if type==0
  %% inputs=[beta,gamma,H], dispersion relation is
  %% F(k)=( beta*k^4+gamma )*k*tanh(k*H)-1=0;
  %% put L=beta^.2, k=K/L
  %%  => ( K^5+(gamma/L)*K )*tanh(K*H/L)-1=0;
  bet=inputs(:,1);
  L=bet.^.2;
  gam_on_L=inputs(:,2)./L;
  HonL=inputs(:,3)./L;
else
  %% inputs=[h,period,H];
  h=inputs(:,1);
  period=inputs(:,2);
  H=inputs(:,3);
  %%
  pram=NDphyspram(0);
  E=pram(1); %Pa
  g=pram(2); %m/s^2
  rho=pram(3); %kg/m^3
  rho_ice=pram(4); %kg/m^3
  nu=pram(5);
  %%
  om=2*pi./period;
  D=E*h.^3/12/(1-nu^2);
  L=( D./rho./om.^2 ).^.2;
  sig=(rho_ice*h)./(rho.*L);
  lam=g./(L.*om.^2);

  %% dispersion relation is
  %% F(k) = ( beta*k^4+gamma )*k*tanh(k*H)-1 = 0,
  %%  where beta = D/rho/om^2,
  %%  gamma = g/om^2 - rho_ice*h/rho;

  gam_on_L=lam-sig;
  HonL=H./L;
end

%% CALCULATE OUTPUTS:
M=length(L);
for j=1:M
  k_cell{jnz(j)}=RTS_ice_roots(gam_on_L(j),HonL(j),N)/L(j);
end
outputs(jnz,:)=[gam_on_L,HonL,L];