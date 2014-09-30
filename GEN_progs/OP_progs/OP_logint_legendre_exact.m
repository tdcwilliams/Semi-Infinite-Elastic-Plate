function Y=OP_logint_legendre_exact(xx,Ngl,do_test)
%% CALL: Y=OP_logint_legendre_exact(xx,Ngl,do_test)
%% xx a column vector;
%% does integral
%%  \int_{-1}^1[ log|xx-x0|*[P_n(x0)/h_n] ]dx0,
%%    where xx\in[-1,1] and n=0..Ngl;
%% P_n are legendre polynomials,
%% h_n=\int_{-1}^1[P_n(t)]^2dt;
%% the answers for each n go into the (n+1)-th column;
%% enter a 3rd argument to do a test.

do_test=(nargin>2);
Y=zeros(length(xx),Ngl+1);
hn=2./(1+2*(0:Ngl));

%% Calculate integral with P_0=1:
n=0;
I0=0*xx;
jp=find(xx~=1);
I0(jp)=(1-xx(jp)).*(log(1-xx(jp))-1);
%%
jm=find(xx~=-1);
I0(jm)=I0(jm)+(1+xx(jm)).*(log(1+xx(jm))-1);
Y(:,n+1)=I0;

if Ngl>0
  %% Calculate integral with P_1=x:
  n=1;
  I1=xx.*I0;
  I1(jp)=I1(jp)+(1-xx(jp)).^2.*(2*log(1-xx(jp))-1)/4;
  I1(jm)=I1(jm)-(1+xx(jm)).^2.*(2*log(1+xx(jm))-1)/4;
  Y(:,n+1)=I1;
  %% Calculate remaining integrals with recurrence relation:
  for n=2:Ngl
    In=(2*n-1)/(n+1)*xx.*I1 + (2-n)/(n+1)*I0 + 2/3*(n==2);
    Y(:,n+1)=In;
    I0=I1;
    I1=In;
  end
end
Y=Y*diag(1./hn);

if do_test%% do a test -try to approx log|x-t| on [-1 1]
          %% with legendre poly's:
  x=xx(1);
  Nterms=Ngl;
  an=OP_logint_legendre_exact(x,Nterms);
  tt=(-1:.001:1)';
  y=log(abs(x-tt));
  plot(tt,y), hold on;
  %%
  y_ap=OP_interp_legendre(tt,an);
  plot(tt,y_ap,'--r'), hold off;
end