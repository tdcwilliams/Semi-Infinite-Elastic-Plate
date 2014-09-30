function [s,ds]=GEN_arc_length(FXN,t0,t,scaler,varargin)
%% CALL: [s,ds]=GEN_arc_length(FXN,t0,t,scaler,varargin)
%% calculates the arc length
%%  of a parametrically-defined curve
%%  as the parameter varies between t0 & t
%% FXN=@f, where [dr,r]=feval(FXN,t,varargin{:})
%%  dr(1,:)=dx, dr(2,:)=dy, r(1,:)=x, r(2,:)=y.

M_sc=diag(scaler);
dr=M_sc*feval(FXN,t,varargin{:});
ds=sqrt( dr(1,:).^2+dr(2,:).^2 );
%%
Ngl=0;%% initialise no of Gauss-legendre points
      %% (gets looped up):
tol=1e-12;%% tolerance in relative error desired
err=2*tol;%% initialise relative error
s0=0;
%%
while err>tol
  Ngl=Ngl+50;
  [tgl,wgl]=OP_numint_legendre(Ngl,[t0 t]);
  dr=M_sc*feval(FXN,tgl,varargin{:});
  ds0=sqrt( dr(1,:).^2+dr(2,:).^2 );
  s=ds0*wgl;
  %%
  err=abs(1-s0/s);
  s0=s;
end