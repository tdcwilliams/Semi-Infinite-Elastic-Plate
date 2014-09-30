function [An,hn,W]=GEN_inprod_exp(tj,wj,Nterms,fj)

%%
DO_TEST=0;
if nargin<2
  Nterms=tj;
  tj=[];
  DO_TEST=1;
elseif nargin<3
  Nterms=round(length(tj)/2)-1;
end
hn=2*ones(2*Nterms+1,1);

if isempty(tj)
  %% get points with usual exponential integration rule,
  %% using 2*(N+1) evenly spaced points:
  [tj,wj]=GEN_numint_exp(Nterms+1);
end

nvec=(-Nterms:Nterms);
%[size(nvec),size((pi*tj).')]
W=exp(-i*pi*tj*nvec);
Exp=W.';

if nargin==4
  An=Exp*(.5*wj.*fj);
else
  An=Exp*diag(wj/2);
end

if DO_TEST
  np=50;
  tt=-1+(0:np)'*2/np;
  %%
  k=4*pi;
  %k=3;
  tau=tt.^2+tt-1;
  tau_j=tj.^2+tj-1;
  ff=3*cos(k*tau);
  fj=3*cos(k*tau_j);
  fn=An*fj;
  fap=GEN_interp_exp(tt,fn);
  plot(tt,[real(ff),imag(ff)]), hold on;
  plot(tt,[real(fap),imag(fap)],'--r'), hold off;
end