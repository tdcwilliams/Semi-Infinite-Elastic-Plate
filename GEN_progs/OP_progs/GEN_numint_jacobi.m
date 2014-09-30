function [x,w]=GEN_numint_jacobi(alf,bet,N)

disp('please use OP_numint_jacobi.m');

MAXIT=10;
EPS=3e-12;

% Initial approximations to the roots:
alfbet=alf+bet;
for j=1:N
z_unf(j,1)=cos( pi*( (N+1-j)-0.25+0.5*alf)/(N+0.5*(alfbet+1)) );
end

z=z_unf;
z1=z;
j_unf=1:N;
p1=0*z;
p2=0*z;
p3=0*z;
pp=0*z;
% Newton’s method carried out simultaneously on the roots.
for its=1:MAXIT
  temp=2+alfbet;
  p1(j_unf)=(alf-bet+temp*z_unf)/2;
  p2(j_unf)=1;

  % Loop up the recurrence relation to get
  % the Jacobi polynomials evaluated at z
  for j=2:N
    temp=2*j+alfbet;%%

    %a=2*j*(j+alfbet)*temp;
    %temp=temp+2;
    a=2*j*(j+alfbet)*(temp-2);%%
    b=(temp-1)*(alf*alf-bet*bet+temp*(temp-2)*z_unf);
    c=2*(j-1+alf)*(j-1+bet)*temp;

    p3(j_unf)=p2(j_unf);
    p2(j_unf)=p1(j_unf);
    p1(j_unf)=(b.*p2(j_unf)-c*p3(j_unf))/a;
  end
  % p1 now contains the desired Legendre polynomials.
  % We next compute pp, the derivatives, by a standard relation
  % involving also p2, the polynomials of one lower order.

  pp(j_unf)=(N*(alf-bet-temp*z_unf).*p1(j_unf)+2*(N+alf)*...
           (N+bet)*p2(j_unf))./(temp*(1-z_unf.^2));

  z1(j_unf)=z_unf;
  z(j_unf)=z_unf-p1(j_unf)./pp(j_unf);
  j_unf=find( (abs(z-z1) > EPS) );
  z_unf=z(j_unf);

  if isempty(j_unf)
    break
  end
end

if (its == MAXIT+1) then
  disp(' ')
  disp('WARNING !! its='),disp(its),
  disp('too many iterations in "GEN_numint_jacobi')
  disp(' ')
end

%OUTPUTS
x=z;
w=gamma(alf+N)*gamma(bet+N)/gamma(N+1)/...
     gamma(N+alf+bet+1)*temp*2^alfbet./(pp.*p2);
