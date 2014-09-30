function [An,hn,Pn_vals]=OP_inprod_legendre(xj,wj,NgP,fj)
%% CALL: [An,hn]=GEN_inprod_legendre(xj,wj,fj)
%% xj\in[-1,1] are the abscissae; wj the weights
%% for Gauss-Legendre integration
%% NB fj is optional - if its not given then An is a matrix
%% such that An*f(xj) are the coefficients
%%
%% A_n = (n+.5)*\int^1_{-1}f(x)P_n(x)dx
%%     ~ (n+.5)*\sum_{j=1}^N [w_j.P_n(x_j).f(x_j)]

if nargin==2
  NgP=length(xj)-1;
end
hn=1./( .5+(0:NgP)' );
Pn_vals=ones(length(xj),NgP+1);

if nargin==4
  P0=0*xj;
  P1=1+P0;
  An=[0;P0];
  An(1)=.5*sum(wj.*fj);

  for n=1:NgP
    Pn=(2-1/n)*xj.*P1-(1-1/n)*P0;
    Pn_vals(:,n+1)=Pn;
    An(n+1)=(n+.5)*sum(wj.*Pn.*fj);
    P0=P1;
    P1=Pn;
  end
else
  %want row vectors because we will calc the rows of a matrix
  xj=xj';
  wj=wj';
  An=zeros(NgP+1,length(xj));
  An(1,:)=.5*wj;

  P0=0*xj;
  P1=1+P0;
  for n=1:NgP
    Pn=(2-1/n)*xj.*P1-(1-1/n)*P0;
    Pn_vals(:,n+1)=Pn';
    An(n+1,:)=(n+.5)*(wj.*Pn);
    P0=P1;
    P1=Pn;
  end
end