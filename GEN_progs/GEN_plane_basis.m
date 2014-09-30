function [v1,v2]  = GEN_plane_basis(nvec)

DO_TEST  = 0;

if ~exist('nvec')
   %nvec  = [0 0 1]'
   %nvec  = [0 1 0]'
   %nvec  = [1 0 0]'
   %nvec  = [1 1 0]'
   %nvec  = [0 1 1]'
   %nvec  = [1 0 1]'
   nvec  = [4 2 1]'
   DO_TEST  = 1;
end

nvec  = nvec/norm(nvec);%%normalise normal vector to plane;
[nmax,jmax] = max(nvec);
%%
j1 = mod(jmax+1,3);
if j1==0
   j1 = 3;
end
j2 = mod(jmax+2,3);
if j2==0
   j2 = 3;
end
n1 = nvec(j1);
n2 = nvec(j2);
%%
v1 = zeros(3,1);
v2 = zeros(3,1);

if n1==0
   v1(j1)   = 1;
   vv       = GEN_plane_basis2d([n2 nmax]);
   v2(j2)   = vv(1);
   v2(jmax) = vv(2);
elseif n2==0
   v1(j2)   = 1;
   vv       = GEN_plane_basis2d([n1 nmax]);
   v2(j1)   = vv(1);
   v2(jmax) = vv(2);
else
   v1(jmax) = 1;
   v1(j1)   = -nmax/n1;
   v2(jmax) = 1;
   v2(j2)   = -nmax/n2;
end

v1 = v1/norm(v1);
v2 = v2-(v2'*v1)*v1;
v2 = v2/norm(v2);

if DO_TEST
   %%check orthogonality;
   v1
   dp1   = nvec'*v1
   %%
   v2
   dp2   = nvec'*v2
   %%
   dp12  = v1'*v2
end
