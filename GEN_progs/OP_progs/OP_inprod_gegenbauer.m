function [An,hn,Cn_vals]=...
            OP_inprod_gegenbauer(xj,wj,alpha,NgP,fj)

%  NgP=length(xj);
%  if length(alpha)==2
%    NgP=alpha(2);
%    alpha(2)=[];
%  end

hn=zeros(NgP+1,1);
Cn_vals=ones(length(xj),NgP+1);

if nargin==5
  An=zeros(NgP+1,1);
  p0=0;
  p1=1;

  beta_n=pi*gamma(2*alpha)*exp( (1-2*alpha)*log(2) ...
              - 2*log(gamma(alpha)) );
  An(1)=alpha*sum(wj.*fj)/beta_n;
  hn(1)=beta_n/alpha;

  for n=1:NgP
    pn=2*(n-1+alpha)/n*xj.*p1...
        -(n+2*alpha-2)/n*p0;
    Cn_vals(:,n+1)=pn;
    beta_n=(n-1+2*alpha)*beta_n/n;
    An(n+1)=(n+alpha)*sum(wj.*pn.*fj)/beta_n;
    hn(n+1)=beta_n/(n+alpha);
    p0=p1;
    p1=pn;
  end
else
  xj=xj.';
  wj=wj.';
  An=zeros(NgP+1,length(xj));

  beta_n=pi*gamma(2*alpha)*exp( (1-2*alpha)*log(2) ...
	  - 2*log(gamma(alpha)) );
  An(1,:)=alpha/beta_n*wj;
  hn(1)=beta_n/alpha;

  p0=0*xj;
  p1=1+p0;
  for n=1:NgP
    pn=2*(n-1+alpha)/n*xj.*p1...
        -(n+2*alpha-2)/n*p0;
    Cn_vals(:,n+1)=pn';
    beta_n=(n-1+2*alpha)*beta_n/n;
    An(n+1,:)=(n+alpha)/beta_n*wj.*pn;
    hn(n+1)=beta_n/(n+alpha);
    p0=p1;
    p1=pn;
  end
end
