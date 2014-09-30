function [f,hn]=OP_interp_gegenbauer(tt,alpha,An)

beta_n=pi*gamma(2*alpha)*exp( (1-2*alpha)*log(2) ...
              - 2*log(gamma(alpha)) );
hn(1)=beta_n/alpha;

if ~iscell(An)
  NgC=length(An)-1;
  f=0*tt;
  if isempty(An)
    return
  end

  C0=f;
  C1=1+f;
  f=f+An(1)*C1;

  for its=1:NgC
    Cn=2*(its+alpha-1)/its*tt.*C1-(its+2*alpha-2)/its*C0;
    f=f+An(its+1)*Cn;
    C0=C1;
    C1=Cn;
    hn(its+1)=beta_n/(its+alpha);
  end
else
  NgC=An{1};
  f=zeros(length(tt),NgC+1);
  if isempty(f)
    return
  end

  C0=zeros(length(tt),1);
  C1=1+C0;
  f(:,1)=C1;

  for its=1:NgC
    Cn=2*(its+alpha-1)/its*tt.*C1-(its+2*alpha-2)/its*C0;
    f(:,its+1)=Cn;
    C0=C1;
    C1=Cn;
    hn(its+1)=beta_n/(its+alpha);
  end
end
