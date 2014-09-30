function Md=OP_diff_chebyshev(Ngc,Lc)

want_coeffs=iscell(Ngc);
if want_coeffs
  An=Ngc{1};
  Ngc=length(An)-1;
end

Md=zeros(Ngc,Ngc);

nn=2*(1:Ngc);

for j=1:2:Ngc
  jj=j:Ngc;
  Md=Md+diag(nn(jj),j-1);
end

Md(1,:)=Md(1,:)/2;
Md=[0*Md(:,1),Md];

if nargin==2
  Md=Md/Lc;
end

if want_coeffs
  Md=Md*An;
end