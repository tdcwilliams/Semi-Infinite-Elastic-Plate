function M=GEN_iprules_ice(gam1,gam2,prams,WANT_EDGE_TERMS)

if nargin==3%%sometimes handy to leave the edge terms off:
  WANT_EDGE_TERMS=1;
end

%%prams={DD,lam,sigsig,H}
DD=prams{1};
lam=prams{2};
sigsig=prams{3};
H=prams{4};
%%
D1=DD(1);
D2=DD(2);
sig1=sigsig(1);
sig2=sigsig(2);
[Gam2,Gam1]=meshgrid(gam2,gam1);
if WANT_EDGE_TERMS
  M=-D1*(Gam1.^2+Gam2.^2);
else
  M=zeros(length(gam1),length(gam2));
end
%%
Lam11=D1*gam1.^4+lam-sig1;
Lam12=D1*Gam2.^4+lam-sig1;
Lam22=D2*Gam2.^4+lam-sig2;
f12=Lam22-Lam12;
S=Gam1+Gam2;
D=Gam1-Gam2;
jp=find(D~=0);
M(jp)=M(jp)+f12(jp)./S(jp)./D(jp);
%%
j0=find(D==0);
if ~isempty(j0)
  hr=D1.^(1/3);
  mu=sig1/hr;
  alpBGzz1=calc_res({D1,lam,sig1,H},gam1).*gam1;
  [Gam2,AlpBGzz1]=meshgrid(gam2,alpBGzz1);
  M(j0)=M(j0)-.5./AlpBGzz1(j0);
end
M=diag(1./Lam11)*M*diag(1./Lam22(1,:));

function y=calc_res(Z2,gamma)
%% y=calc_res(Z2,hr,gamma)=Res(1/f(K),gamma_n),
%% where gamma_n is a root of the dispersion relation
%% f=1/K/tanh(KH)-(Dr*K^4+lam-mr*mu);
%% Z2={lam,mu,H}.
Dr=Z2{1};
lam=Z2{2};
sig=Z2{3};
H=Z2{4};
%%
Gam=Dr*gamma.^4+lam-sig;
Gampr=Gam+4*Dr*gamma.^4;
denom=H*(Gam.^2.*gamma.^2-1)+Gampr;
y=-gamma./denom;
