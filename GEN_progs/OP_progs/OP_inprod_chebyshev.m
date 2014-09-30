function [An,hn,Tn_vals]   = OP_inprod_chebyshev(xj,wj,NgP,fj)
%% CALL: [An,hn,Tn_vals]   = OP_inprod_chebyshev(xj,wj,NgP,fj)

if ~exist('NgP')
   NgP = length(xj)-1;
else
   if isempty(NgP)
      NgP = length(xj)-1;
   end
end

hn       = pi/2*ones(NgP+1,1);
hn(1)    = pi;
Tn_vals  = ones(length(xj),NgP+1);

%nvec=(1:NgP)';
%xj=cos(pi/NgP*(nvec-.5));cos(NgP*acos(xj))
%wj=pi/NgP;

if nargin==4
  An  = zeros(NgP+1,1);
  p0  = 0;
  p1  = 1;

  An(1)  = sum(wj.*fj)/hn(1);

  for n=1:NgP
    pn               = ( 2-(n==1) )*xj.*p1 - p0;%%T_1(x)=x, not 2x
    Tn_vals(:,n+1)   = pn;
    An(n+1)          = sum(wj.*pn.*fj)/hn(n+1);
    p0               = p1;
    p1               = pn;
  end
else
  xj  = xj.';
  wj  = wj.';
  An  = zeros(NgP+1,length(xj));

  An(1,:)   = wj/hn(1);

  p0  = 0*xj;
  p1  = 1+p0;
  for n=1:NgP
    pn               = ( 2-(n==1) )*xj.*p1-p0;
    Tn_vals(:,n+1)   = pn';
    An(n+1,:)        = wj.*pn/hn(n+1);
    p0               = p1;
    p1               = pn;
  end
end
