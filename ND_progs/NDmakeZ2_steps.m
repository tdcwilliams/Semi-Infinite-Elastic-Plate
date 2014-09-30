function [Rts,rts,HH,el,lam,sigsig,L,rts_approx]=NDmakeZ2_steps(Z,hh,NN)
%% Z={T,theta_inc,HH_dim}


%% get impt parameters:
[Rts,HH,th,lam,sigsig,L]=ND2_steps(Z,hh,NN);

%[Rts{1}(1),Rts{2}(1)],pause

%% sort out roots:
gam0=Rts{1}(1);
el=gam0*sin(th);
for j=1:2
  gam=Rts{j};
  alp=sqrt(gam.^2-el^2);
  rts{j}=sign((imag(alp)>=0)-.5).*alp;
  %rts_approx{j}=( 1:NNN(j) )'*pi*i/HHsig(j);
end

%  if HHsig(3)==0
%    rts_approx{3}=[];
%  end