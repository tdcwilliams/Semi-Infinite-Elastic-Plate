function y=OP_interp_jacobi(z,alf,bet,An)

alfbet=alf+bet;

if ~iscell(An)
  N=length(An)-1;
  y=0*z;
  if N<0
    return;
  end

  temp=2+alfbet;
  p2=1;
  y=y+An(1)*p2;

  if N>0
    p1=(alf-bet+temp*z)/2;
    y=y+p1*An(2);
  end

  % Loop up the recurrence relation to get
  % the Jacobi polynomials evaluated at z
  for j=2:N
    temp=2*j+alfbet;%%

    %a=2*j*(j+alfbet)*temp;
    %temp=temp+2;
    a=2*j*(j+alfbet)*(temp-2);%%
    b=(temp-1)*(alf*alf-bet*bet+temp*(temp-2)*z);
    c=2*(j-1+alf)*(j-1+bet)*temp;

    p3=p2;
    p2=p1;
    p1=(b.*p2-c*p3)/a;
    y=y+An(j+1)*p1;
  end
else
  N=An{1};
  y=zeros(length(z),N+1);
  if N<0
    return;
  end

  temp=2+alfbet;
  p2=1+0*z;
  y(:,1)=p2;

  if N>0
    p1=(alf-bet+temp*z)/2;
    y(:,2)=p1;
  end

  % Loop up the recurrence relation to get
  % the Jacobi polynomials evaluated at z
  for j=2:N
    temp=2*j+alfbet;%%

    a=2*j*(j+alfbet)*(temp-2);%%
    b=(temp-1)*(alf*alf-bet*bet+temp*(temp-2)*z);
    c=2*(j-1+alf)*(j-1+bet)*temp;

    p3=p2;
    p2=p1;
    p1=(b.*p2-c*p3)/a;
    y(:,j+1)=p1;
  end

end