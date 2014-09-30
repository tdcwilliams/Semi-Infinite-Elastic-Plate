function Y=GEN_intmat_exp_legendre(kk,na)
%% CALL: Y=GEN_intmat_exp_legendre(k,[Ngl,a])
%% does integral
%% \int_0^a[ exp(i*k*x)*P_n(2*x-1)dx
%% NB only stable for Arg[k]\in[-pi,pi]

Ngl=na(1);
a=na(2);

Y=zeros(length(kk),Ngl+1);
% Calculate integral with P_0=1:
n=0;
I0=0*kk+a;%?+2;
jp=find(kk);
k=a*kk(jp)/2;
%I0(jp)=2*sin(k)./k;
I0(jp)=a/2i*(exp(2i*k)-1)./k;
Y(:,n+1)=I0;
I0=I0(jp);

if Ngl>0 & ~isempty(jp)
%% NB don't need integrals with P_n (n>0) if k==0
%% because of orthogonality relation.

  %% Calculate integral with P_1=x:
  n=1;
  %I1=-i./k.^2.*(k.*cos(k)-sin(k));
  I1=-a/2./k.^2.*((i*k-1).*exp(2i*k)+i*k+1);
  Y(jp,n+1)=I1;

  %% Calculate remaining integrals with recurrence relation:
  for n=2:Ngl
    In=(2*n-1)*i*I1./k+I0;
    Y(jp,n+1)=In;
    I0=I1;
    I1=In;
  end
end

if 0%%test results numerically using quad:
  Mt=zeros(length(kk),Ngl+1);
  for j=1:length(kk)
    for n=0:Ngl
      Mt(j,n+1)=quad(@testig,0,a,[],[],kk(j),[n a]);
    end
  end
  Mt
end

function y=testig(x,k,na)

n=na(1);
a=na(2);
An=zeros(n+1,1);
An(n+1)=1;
Pn=GEN_interp_legendre(2*x/a-1,An);
y=Pn.*exp(i*k*x);