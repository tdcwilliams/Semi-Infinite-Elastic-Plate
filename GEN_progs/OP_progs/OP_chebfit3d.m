function [Fcoeffs,Lims]=OP_chebfit3d(x,y,z,Fvals,NNN,DO_TEST)

if nargin==5
   DO_TEST  = 0;
end

if isempty(x)
   DO_TEST  = 1;
   nx       = 19;
   x        = (0:nx)'/nx*2;
   ny       = 22;
   y        = -1+(0:ny)'/ny*4;
   nz       = 32;
   z        = -1+(0:nz)'/nz*4;
   for j=1:nz+1
      z0             = z(j);
      if 0
         Fvals(:,:,j)   = (x.^3-1+x)*(y.^4+3*y-1.5)'*...
                            (z0.^3+1-z0);
      else
         Fvals(:,:,j)   = sin(x)*cos(2*y)'*...
                            cos(z0.^2+1);
      end
   end
end

Lims  = [min(x),max(x),min(y),max(y),min(z),max(z)];
nx    = length(x)-1;
ny    = length(y)-1;
nz    = length(z)-1;
nnn   = [nx,ny,nz];
for j=1:3
   NNN(j)   = min(nnn(j),NNN(j));
end
%%
for j=1:length(z)
   F0(:,:,j)   = OP_chebfit2d( x,y,Fvals(:,:,j),NNN(1:2) );
end
F0a   = permute(F0,[3 1 2]);%% reverse is [2 3 1]
%%
for j=1:NNN(2)+1
   F0b(:,:,j)  = OP_chebfit( NNN(3),z,F0a(:,:,j) );
end

Fcoeffs  = permute(F0b,[2 3 1]);

if DO_TEST
   nt = 2;
   if nt==1%% TEST & LOOK AT COEFFS (plot vs x):
      np    = 100; 
      xt   = Lims(1)+(0:np)'/np*(Lims(2)-Lims(1));
      fap1  = OP_chebinterp3d(xt,y,z,Fcoeffs,Lims);
      %%
      for j=1:ny+1
         disp(['y = ',num2str(y(j))]);
         for r=1:nz+1
            disp(['z = ',num2str(z(r))]);
            plot( x,squeeze(Fvals(:,j,r)),'xk' ), hold on;
            %%
            plot( xt,squeeze(fap1(:,j,r)),'r' ), hold off;
            pause;
         end
         disp(' ');
      end
   elseif nt==2%% TEST & LOOK AT COEFFS (plot vs y): 
      np    = 100; 
      yt   = Lims(3)+(0:np)'/np*(Lims(4)-Lims(3));
      fap1  = OP_chebinterp3d(x,yt,z,Fcoeffs,Lims);
      %%
      for j=1:nx+1
         disp(['x = ',num2str(x(j))]);
         for r=1:nz+1
            disp(['z = ',num2str(z(r))]);
            plot( y,squeeze(Fvals(j,:,r)),'xk' ), hold on;
            %%
            plot( yt,squeeze(fap1(j,:,r)),'r' ), hold off;
            pause;
         end
         disp(' ');
      end
   else%% TEST & LOOK AT COEFFS (plot vs z):
      np    = 100; 
      zt   = Lims(5)+(0:np)'/np*(Lims(6)-Lims(5));
      fap1  = OP_chebinterp3d(x,y,zt,Fcoeffs,Lims);
      %%
      for j=1:nx+1
         disp(['x = ',num2str(x(j))]);
         for r=1:ny+1
            disp(['y = ',num2str(y(r))]);
            plot( z,squeeze(Fvals(j,r,:)),'xk' ), hold on;
            %%
            plot( zt,squeeze(fap1(j,r,:)),'r' ), hold off;
            pause;
         end
         disp(' ');
      end
   end
end
