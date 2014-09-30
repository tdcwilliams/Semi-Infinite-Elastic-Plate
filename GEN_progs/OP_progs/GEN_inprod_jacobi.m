function [An,hn]=GEN_inprod_jacobi(fj,xj,wj,alpha,beta)

disp('please use OP_inprod_jacobi.m');

NgP=length(fj);
p0=1;
beta_n=exp( (1+alpha+beta)*log(2) )*gamma(1+alpha)...
          *gamma(1+beta)/gamma(1+alpha+beta);
An(1)=(alpha+beta+1)*sum(wj.*fj)/beta_n;
hn=beta_n/(alpha+beta+1)+zeros(NgP+1,1);

n=1;
p1=(alpha-beta+(alpha+beta+2)*xj)/2
beta_n=(n+alpha)*(n+beta)/n/(n+alpha+beta)*beta_n;
An(2)=(2*n+alpha+beta+1)*sum(wj.*p1.*fj)/beta_n;
hn(n+1)=beta_n/(2*n+alpha+beta+1);

for n=2:NgP
  r1n=(n+alpha+beta);%2b continued
  r2n=(r1n+n-1)*(alpha^2-beta^2);
  r3n=(r1n+n-2)*(r1n+n-1)*(r1n+n);
  r4n=2*(n+alpha-1)*(n+beta-1)*(n+r1n);
  r1n=2*n*r1n*(r1n+n-2);

  pn=( (r2n+r3n*xj).*p1 -r4n*p0 )/r1n;
  beta_n=(n+alpha)*(n+beta)/n/(n+alpha+beta)*beta_n;
  An(n+1)=(2*n+alpha+beta+1)*sum(wj.*pn.*fj)/beta_n;
  hn(n+1)=beta_n/(2*n+alpha+beta+1);
  p0=p1;
  p1=pn;
end
