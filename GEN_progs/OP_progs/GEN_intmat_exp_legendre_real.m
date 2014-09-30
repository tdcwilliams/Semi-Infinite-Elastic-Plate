function [Y,Mt]=GEN_intmat_exp_legendre_real(kk,ab,Ngl)
%% CALL: Y=GEN_intmat_exp_legendre_real(k,[a b],Ngl)
%% does integral
%% \int_a^b[ exp(i*k*x)*P_n(x)dx
%% where general formula is only stable for real k.

do_test=(nargout==2);

%% allow for arbitrary integration range [a,b]
a=ab(1);
b=ab(2);
kfac=(b-a)/2;
Y=zeros(length(kk),Ngl+1);
kc=1.1;

%% DO k=0 EXACTLY
j0=find(kk==0);
Y(j0,1)=b-a;

%% DO SMALL k's (~=0) NUMERICALLY
%% AS PROBABLY NOT OSCILLATING MUCH
%% AND EXACT FORMULA CAN BE UNSTABLE FOR SMALL k~=0 & LARGE n:
jm=find(abs(kk)<=kc & kk~=0);
[xgl,wgl]=GEN_numint_legendre(a,b,Ngl);
tgl=-1+2*(xgl-a)/(b-a);
E=exp(i*kk(jm)*xgl');
M=GEN_AD(E/i,i*wgl);%%=E*diag(wgl);
Y(jm,:)=M*GEN_interp_legendre(tgl,{Ngl});

%% DO LARGER k's EXACTLY:
jp=find(abs(kk)>kc);
fac=kfac*exp(i*kk(jp)*(a+b)/2);
%% use relation I_n(a,b,k)=fac*I_n(-1,1,kfac*k)

% Calculate integral with P_0=1:
n=0;
%I0=0*kk+b-a;
k=kfac*kk(jp);
I0=2*fac.*sin(k)./k;
Y(jp,n+1)=I0;
%I0=I0(jp);

if Ngl>0 & ~isempty(jp)
%% NB don't need integrals with P_n (n>0) if k==0
%% because of orthogonality relation.

  %% Calculate integral with P_1=x:
  n=1;
  I1=-2i*fac./k.^2.*(k.*cos(k)-sin(k));
  Y(jp,n+1)=I1;

  %% Calculate remaining integrals with recurrence relation:
  for n=2:Ngl
    In=(2*n-1)*i*I1./k+I0;
    Y(jp,n+1)=In;
    I0=I1;
    I1=In;
  end
end

if do_test%%test results numerically using quad:
  Mt=0*k;
  for j=1:length(kk)
    Mt(j)=quad(@testig,a,b,[],[],kk(j),ab,Ngl);
  end
elseif 0
  Mt=zeros(length(kk),Ngl+1);
  for j=1:length(kk)
    for n=0:Ngl
      Mt(j,n+1)=quad(@testig,a,b,[],[],kk(j),ab,n);
    end
  end
end

function y=testig(x,k,ab,n)

a=ab(1);
b=ab(2);
t=-1+2*(x-a)/(b-a);

An=zeros(n+1,1);
An(n+1)=1;
Pn=GEN_interp_legendre(t,An);
y=Pn.*exp(i*k*x);