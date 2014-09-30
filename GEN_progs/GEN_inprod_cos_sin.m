function [An,hn,W]=GEN_inprod_cos_sin(tj,wj,Nterms,fj)

hn=ones(2*Nterms+1,1);
hn(1)=2;
%%
DO_TEST=0;
if nargin==1
  Nterms=tj;
  tj=[];
  DO_TEST=1;
end

if isempty(tj)%%get points with Gauss-Legendre:
  ppw=4;
  Nint=ppw*Nterms;
  [tj,wj]=OP_numint_legendre(Nint);
  if DO_TEST
    fj=tj;
  end
end

W=ones(length(tj),2*Nterms+1);
cos_th_j=cos(pi*tj);
sin_th_j=sin(pi*tj);

if nargin==4
  An_cos=zeros(Nterms+1,1);
  Bn_sin=zeros(Nterms,1);
  p0=0;
  p1=1;
  q1=0;
  %%
  An_cos(1)=sum(wj.*fj)/2;
  for n=1:Nterms
    pn=( 2-(n==1) )*cos_th_j.*p1 - p0;
         %% n=1->cos(th), else
         %% cos(n*th)=2*cos(th)*cos((n-1)*th)-cos((n-2)*th)
    An_cos(n+1)=sum(wj.*pn.*fj);
    W(:,n+1)=pn;
    %%
    qn=cos_th_j.*q1+sin_th_j.*p1;
         %% sin(n*th)=cos(th)*sin((n-1)*th)
         %%               +sin(th)*cos((n-1)*th)
    Bn_sin(n)=sum(wj.*qn.*fj);
    W(:,n+1+Nterms)=qn;
    %%
    p0=p1;
    p1=pn;
    q1=qn;
  end
else
  cos_th_j=cos_th_j.';
  wj=wj.';
  An_cos=zeros(Nterms+1,length(cos_th_j));
  %%
  sin_th_j=sin_th_j.';
  Bn_sin=zeros(Nterms,length(cos_th_j));
  %%
  An_cos(1,:)=wj/2;
  p0=0*cos_th_j;
  p1=1+p0;
  q1=p0;
  %%
  for n=1:Nterms
    pn=( 2-(n==1) )*cos_th_j.*p1-p0;
    An_cos(n+1,:)=wj.*pn;
    W(:,n+1)=pn;
    %%
    qn=cos_th_j.*q1+sin_th_j.*p1;
    Bn_sin(n,:)=wj.*qn;
    W(:,n+1+Nterms)=qn;
    %%
    p0=p1;
    p1=pn;
    q1=qn;
  end
end

An=[An_cos;Bn_sin];
if DO_TEST
  fap=GEN_interp_cos_sin(tj,An*fj);
  plot(tj,fj), hold on;
  plot(tj,fap,'--r'), hold off;
end