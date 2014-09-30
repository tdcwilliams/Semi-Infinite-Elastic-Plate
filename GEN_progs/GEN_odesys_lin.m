function an=GEN_odesys_lin(Af_fxn,Nterms,interval,varargin)
%% CALL: an=GEN_odesys_lin(Af_fxn,Nterms,interval,varargin)
%%  finds general solution to Y'-A*Y=f;
%%   gives solution as
%%    Y_j(y)=\sum_{n=0}^Nterms a_n^(j)T_n(t),
%%     where the a_0^(j) are arbitrary constants,
%%      t=-1+2*(y-a)/(b-a) => y=a+(b-a)*(t+1)/2;
%% interval =[a,b] - domain eqn is over;
%% include extra arguments for Af_fxn after interval

DO_TEST=0;
if nargin==0%% USE SOME TEST INPUTS:
  Nterms=5;
%    Af_fxn=@test_fn1;
  Af_fxn=@test_fn2;
  DO_TEST=1;
end

if nargin<3
  a=-1;
  b=1;
else
  a=interval(1);
  b=interval(2);
end

Nint=max(50,2*Nterms);
alpU=1;
[tt,ww]=OP_numint_gegenbauer(alpU,Nint);
yy=a+.5*(b-a)*(1+tt);
%%
[ipU,hn]=OP_inprod_gegenbauer(tt,ww,alpU,Nterms-1);
Tn_vals=OP_interp_chebyshev(tt,{Nterms});

[Amatrix,forcing]=feval(Af_fxn,yy,varargin{:});
M=round(sqrt( size(Amatrix,2) ));
%%
nvec=(1:Nterms)';
ipU2=diag(.5*(b-a)./nvec)*ipU;
%%
KMAT=zeros(M*Nterms,M*(Nterms+1));
FORCING=KMAT(:,1);
%%
for j=1:M
  j_rows=nvec+(j-1)*Nterms;
  FORCING(j_rows)=ipU2*forcing(:,j);
  %%
  for r=1:M
    j_cols=[1;nvec+1]+(r-1)*(Nterms+1);
    s=j+(r-1)*M;
    KMAT(j_rows,j_cols)=ipU2*diag(Amatrix(:,s))*Tn_vals;
  end
end

j0=1:Nterms+1:M*(Nterms+1);
FORCING=[FORCING,KMAT(:,j0)];
KMAT(:,j0)=[];
an0=( eye(M*Nterms)-KMAT )\FORCING;

II=[zeros(M,1),eye(M)];
an=zeros(size([II;an0]));
for j=1:M
  s=j0(j);
  an(s,:)=II(j,:);
  an(s+nvec,:)=an0(nvec+(j-1)*Nterms,:);
end

if DO_TEST
  for j=1:M
    jj=(1:Nterms+1)+(j-1)*(Nterms+1);
    subplot(1,M,j), plot(yy,Tn_vals*an(jj,:));
  end
  c1=Tn_vals(1,:)*an(jj,M+1);
  c0=exp(sin(yy(1)));
  hold on, plot(yy,c1/c0*exp(sin(yy)),'--k');
  plot(yy,sin(yy),'--g');
  %%
  c1=Tn_vals(1,:)*an(1:Nterms+1,2);
  c0=exp(-cos(yy(1)));
  subplot(1,M,1), hold on;
  plot(yy,c1/c0*exp(-cos(yy)),'--r');
  plot(yy,sin(yy),'--g');
end

function [A,f]=test_fn1(y,varargin)

A=cos(y);
f=0*y;
%% general solution is c*e^{sin(x)}

function [A,f]=test_fn2(y,varargin)

z=0*y;
f=[cos(y)-sin(y).^2, cos(y)-cos(y).*sin(y)];
A=[sin(y),z,z,cos(y)];

%% general solution is [c1*e^{-cos(x)};c2*e^{sin(x)}]