function [DET,dDET]=GEN_diff_det33(M,dM);
%% differentiates the determinant of a 3x3 matrix,
%%   given the matrix and its derivative;

DO_TEST=0;
if nargin==0
  kk=(0:.1:1)';
  mat=[1 2 3;4 5 6;1 8 9];
  for j=1:3
    for r=1:3
      M{j,r}=mat(j,r)*kk.^2;
      dM{j,r}=2*mat(j,r)*kk;
    end
  end
  DO_TEST=1;
end

a=M{1,1};
b=M{1,2};
c=M{1,3};
d=M{2,1};
f=M{2,2};
g=M{2,3};
h=M{3,1};
j=M{3,2};
k=M{3,3};
%%
da=dM{1,1};
db=dM{1,2};
dc=dM{1,3};
dd=dM{2,1};
df=dM{2,2};
dg=dM{2,3};
dh=dM{3,1};
dj=dM{3,2};
dk=dM{3,3};
%%
DET = + a.*(f.*k-g.*j) + ...
      - b.*(d.*k-g.*h) + ...
      + c.*(d.*j-f.*h);
%%
dDET = + da.*(f.*k-g.*j) + ...
       - db.*(d.*k-g.*h) + ...
       + dc.*(d.*j-f.*h) + ...
       + a.*(df.*k+k.*df-dg.*j-g.*dj) + ...
       - b.*(dd.*k+d.*dk-dg.*h-g.*dh) + ...
       + c.*(dd.*j+dj.*d-df.*h-f.*dh);

if DO_TEST
  [DET,18*kk.^6]
  [dDET,6*18*kk.^5]
end