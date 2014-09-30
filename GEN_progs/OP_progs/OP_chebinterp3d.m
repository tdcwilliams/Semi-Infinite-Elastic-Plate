function F=OP_chebinterp3d(xx,yy,zz,F_coeffs,Lims);

x_min    = Lims(1);
x_max    = Lims(2);
y_min    = Lims(3);
y_max    = Lims(4);
z_min    = Lims(5);
z_max    = Lims(6);
%%
Ncheb1   = size(F_coeffs,1)-1;
Ncheb2   = size(F_coeffs,2)-1;
Ncheb3   = size(F_coeffs,3)-1;
%%
dx       = x_max - x_min;
tt1      = -1 + 2*(xx-x_min)/dx;
TnVals1  = OP_interp_chebyshev(tt1,{Ncheb1});
%%
dy       = y_max - y_min; 
tt2      = -1 + 2*(yy-y_min)/dy;
TnVals2  = OP_interp_chebyshev(tt2,{Ncheb2 });
%%
dz       = z_max - z_min;
tt3      = -1 + 2*(zz-z_min)/dz;
TnVals3  = OP_interp_chebyshev(tt3,{Ncheb3});
%%
n1 = length(xx);
n2 = length(yy);
n3 = length(zz);
%%
for j=1:Ncheb3+1
   F0(:,:,j)   = TnVals1*F_coeffs(:,:,j)*TnVals2';
end
F0a   = permute(F0,[3 1 2]);
%% n1 x n2 x Ncheb3 -> Ncheb3 x n1 x n2;

for j=1:n2
   F0b(:,:,j)  = TnVals3*F0a(:,:,j);
end
F  = permute(F0b,[2 3 1]);
