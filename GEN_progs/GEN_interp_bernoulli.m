function y=GEN_interp_Bernoulli(x,ck)

do_test=0;
if nargin==0
  do_test=1;
  N=3;
  if 0
    ck=zeros(N+1,1);
    ck(N+1)=1;
    ckt=1;
  else
    ck={N};
    ckt=zeros(N+1,1);
    ckt(N+1)=1;
  end
  np=100;
  x=(0:np)'/np;
end

if ~iscell(ck)
  N=length(ck)-1;
  an0=GEN_bernoulli_nos(N);
  y=0*x+sum(ck.*an0);
  jj=(1:N).';
  ck(1)=[];
  xx=1+0*x;
  for r=1:N
    xx=x.*xx;
    an=jj.*an0(jj+1-r)/r;
    y=y+sum(ck.*an)*xx;
    an0=an;
    ck(1)=[];
    jj(1)=[];
  end
else
  N=ck{1};
  an0=GEN_bernoulli_nos(N);
  y=zeros(length(x),N+1);
  for r=0:N
    y(:,r+1)=an0(r+1);
  end
  %%
  jj=(1:N).';
  xx=1+0*x;
  for s=1:N
    xx=x.*xx;
    an=jj.*an0(jj+1-s)/s;
    for r=s:N
      y(:,r+1)=y(:,r+1)+an(r+1-s)*xx;
    end
    jj(1)=[];
    an0=an;
  end
end

if do_test
  plot(x,y*ckt)
  if N==0
    yt=1+0*x;
    hold on, plot(x,yt,'--r'), hold off;
  elseif N==1
    yt=x-1/2;
    hold on, plot(x,yt,'--r'), hold off;
  elseif N==2
    yt=x.^2-x+1/6;
    hold on, plot(x,yt,'--r'), hold off;
  elseif N==3
    yt=x.^3-3/2*x.^2+1/2*x;
    hold on, plot(x,yt,'--r'), hold off;
  elseif N==4
    yt=x.^4-2*x.^3+x.^2-1/30;
    hold on, plot(x,yt,'--r'), hold off;
  elseif N==5
    yt=x.^5-5/2*x.^4+5/3*x.^3-1/6*x;
    hold on, plot(x,yt,'--r'), hold off;
  elseif N==6
    yt=x.^6-3*x.^5+5/2*x.^4-1/2*x.^2+1/42;
    hold on, plot(x,yt,'--r'), hold off;
  end
end