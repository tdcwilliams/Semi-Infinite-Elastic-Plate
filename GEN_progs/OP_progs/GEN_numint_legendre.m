% this module provides the subroutine GEN_gauleg, which computes
% N-vectors x & w of abscissae and weights
% to use for Gauss-Legendre integration
% given N and the endpoints x1 & x2.
% CALL [x,w]=GEN_gauleg(x1,x2,N)function

function [x,w]=GEN_numint_legendre(x1,x2,N)

disp('please use OP_numint_Legendre.m');

EPS=3e-14;
MAXIT=10;

M=floor((N+1)/2); % The roots are symmetric in the interval,
          % so we only have to ﬁnd half of them.
xm=0.5*(x2+x1);
xl=0.5*(x2-x1);

% Initial approximations to the roots:
zz=cos(pi*((1:M)'-0.25)/(N+0.5));
ppvec=0*zz;
zzvec=ppvec;
x=zeros(N,1);
w=x;

for jj=1:M
  unfinished=1;
  z=zz(jj);
  % Newton’s method carried out individually on the roots.
  for its=1:MAXIT
    p1=1;
    p2=0;
    % Loop up the recurrence relation to get
    % the Legendre polynomials evaluated at z
    if unfinished
      for j=1:N
         p3=p2;
         p2=p1;
         p1=((2*j-1)*z*p2-(j-1)*p3)/j;
      end
      % p1 now contains the desired Legendre polynomials.
      % We next compute pp, the derivatives, by a standard relation
      % involving also p2, the polynomials of one lower order.
      pp=N*(z*p1-p2)/(z*z-1);
      z1=z;
      z=z1-p1/pp;
      unfinished=(abs(z-z1) > EPS);
    else
      zzvec(jj)=z;
      ppvec(jj)=pp;
      break
    end
  end
  if (its == MAXIT+1)
    disp('too many iterations in GEN_gauleg')
  end
end

%OUTPUTS
x(1:M)=xm-xl*zzvec;        % Scale the root to the desired interval
x(N:-1:N-M+1)=xm+xl*zzvec; % Put in its symmetric counterpart.
w(1:M)=2*xl./((1-zzvec.^2).*ppvec.^2);  % Compute the weight
w(N:-1:N-M+1)=w(1:M);                   % and its symmetric counterpart
