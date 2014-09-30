function Y=GEN_interp_cos_sin(tj,An)
%% CALL Y=GEN_interp_cos_sin(tj,An)
%%
%% usually An=[A_n;B_n] =>
%% Y=\sum_{n=0}^N A_n*cos(n*pi*tj)+\sum_{n=1}^N B_n*sin(n*pi*tj);
%%
%% otherwise An=cell{Nterms} =>
%% Y is a matrix s.t. Y*An2 evaluates the series.

DO_TEST=0;
if DO_TEST
  ntst=13;
  ppw=50;
  np=ppw*ntst;
  tj=-1+(0:np-1)'*2/(np-1);
  SIN=1*(ntst>0);
  An_tst=zeros(2*ntst+1,1);
  An_tst(ntst+1+SIN*ntst)=1;
  TESTMAT=1;
  if TESTMAT
    An={ntst};
  else
    An=An_tst;
  end
end

cos_th_j=cos(pi*tj);
sin_th_j=sin(pi*tj);

if iscell(An)==0
  Y=0*tj;
  Nterms=(length(An)-1)/2;
  An_cos=An(1:Nterms+1);
  Bn_sin=An(Nterms+2:2*Nterms+1);
  p0=0;
  p1=1;
  q1=0;
  %%
  Y=An_cos(1);
  for n=1:Nterms
    pn=( 2-(n==1) )*cos_th_j.*p1 - p0;%%T_1(x)=x, not 2x
    Y=Y+An_cos(n+1)*pn;
    %%
    qn=cos_th_j.*q1+sin_th_j.*p1;
    Y=Y+Bn_sin(n)*qn;
    %%
    p0=p1;
    p1=pn;
    q1=qn;
  end
else
  Nterms=An{1};
  Ycos=ones(length(cos_th_j),Nterms+1);
  Ysin=ones(length(cos_th_j),Nterms);
  %%
  p0=0*cos_th_j;
  p1=1+p0;
  q1=p0;
  %%
  for n=1:Nterms
    pn=( 2-(n==1) )*cos_th_j.*p1-p0;
    Ycos(:,n+1)=pn;
    %%
    qn=cos_th_j.*q1+sin_th_j.*p1;
    Ysin(:,n)=qn;
    %%
    p0=p1;
    p1=pn;
    q1=qn;
  end
  Y=[Ycos,Ysin];
end

if DO_TEST
  if SIN
    Ytst=sin(ntst*pi*tj);
  else
    Ytst=cos(ntst*pi*tj);
  end
  if TESTMAT
    Y2=Y*An_tst;
    plot(tj,Y2), hold on;
    plot(tj,Ytst,'--r'), hold off;
  else
    plot(tj,Y), hold on;
    plot(tj,Ytst,'--r'), hold off;
  end
end