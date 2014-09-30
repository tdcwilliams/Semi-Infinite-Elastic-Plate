function Y=GEN_int_log_legendre(xx,si,An)
%% CALL: Y=int_logtimesPnseries(xx,An)
%% xx\in[-1,1]
%% does integral \int_{-1}^1[ log|xx-s|.\sum_{n=0}^Ngl[A_n.P_n(s)] ]ds
Ngl=length(An)-1;
Y=0*xx;
% Calculate integral with P_0=1:
n=0;
I0=Y;
jp=find(xx~=1);
I0(jp)=(1-xx(jp)).*(log(1-xx(jp))-1);
%%
jm=find(xx~=-1);
I0(jm)=I0(jm)+(1+xx(jm)).*(log(1+xx(jm))-1);
Y=Y+An(n+1)*I0;

if Ngl>0
  % Calculate integral with P_1=x:
  n=1;
  I1=xx.*I0;
  I1(jp)=I1(jp)+(1-xx(jp)).^2.*(2*log(1-xx(jp))-1)/4;
  I1(jm)=I1(jm)-(1+xx(jm)).^2.*(2*log(1+xx(jm))-1)/4;
  Y=Y+An(n+1)*I1;
  % Calculate remaining integrals with recurrence relation:
  for n=2:Ngl
    In=(2*n-1)/(n+1)*xx.*I1 + (2-n)/(n+1)*I0 + 2/3*(n==2);
    Y=Y+An(n+1)*In;
    I0=I1;
    I1=In;
  end
end

% allow for si~=1:
Y=si*(Y+2*An(1)*log(si));
