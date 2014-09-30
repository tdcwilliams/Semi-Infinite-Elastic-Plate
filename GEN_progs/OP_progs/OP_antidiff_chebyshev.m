function Mad=OP_antidiff_chebyshev(Ngl,Lc)

Mad=zeros(Ngl+2,Ngl+1);

n=1:Ngl+1;
Mad(2:Ngl+2,:)=diag(1./(2*n));
Mad(2,1)=1;

n=1:Ngl-1;
Mad(2:Ngl,3:Ngl+1)=Mad(2:Ngl,3:Ngl+1)-diag(1./(2*n));

if nargin==2
  Mad=Mad*Lc;
end