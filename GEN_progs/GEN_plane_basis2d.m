function v1 = GEN_plane_basis2d(nvec)

DO_TEST  = 0;

if ~exist('nvec')
   nvec  = [0 1]'
   %nvec  = [1 0]'
   %nvec  = [1 1]'
   DO_TEST  = 1;
end

nvec        = nvec/norm(nvec);%%normalise normal vector to plane;
[nmax,jmax] = max(nvec);
%%
j1 = mod(jmax+1,2);
if j1==0
   j1 = 2;
end
n1 = nvec(j1);
%%
v1 = zeros(2,1);

if n1==0
   v1(j1)   = 1;
else
   v1(jmax) = 1;
   v1(j1)   = -nmax/n1;
end
v1 = v1/norm(v1);

if DO_TEST
   v1
   dp1   = nvec'*v1
end
