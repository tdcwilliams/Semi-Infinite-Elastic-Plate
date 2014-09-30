function [Fcoeffs,Lims]=OP_chebfit2d(x,y,Fvals,NN,DOTEST)

if nargin==4
   DO_TEST  = 0;
end

if isempty(x)
   DO_TEST  = 1;
   nx       = 19;
   x        = (0:nx)'/nx*2;
   ny       = 22;
   y        = -1+(0:ny)'/ny*4;
   if 0
      Fvals = (x.^3-1+x)*(y.^4+3*y-1.5)';
   else
      Fvals = sin(2*pi*x)*cos(y.^2-1)';
   end
end

Lims     = [min(x),max(x),min(y),max(y)];
nx    = length(x)-1;
ny    = length(y)-1;
nn    = [nx,ny];
for j=1:2
   NN(j)   = min(nn(j),NN(j));
end
%%
F0       = OP_chebfit(NN(1),x,Fvals);
F0a      = transpose(F0);
%%
F0b      = OP_chebfit(NN(2),y,F0a);
Fcoeffs  = transpose( F0b );

if DO_TEST
   if 0%% TEST & LOOK AT COEFFS: 
      np    = 100; 
      xt   = Lims(1)+(0:np)'/np*(Lims(2)-Lims(1));
      fap1  = OP_chebinterp2d(xt,y,Fcoeffs,Lims);
      %%
      for j=1:ny+1
         disp(['y = ',num2str(y(j))]);
         %%
         plot( x,Fvals(:,j),'xk' ), hold on;
         %%
         plot( xt,fap1(:,j),'r' ), hold off;
         pause;
      end
   elseif 1%% TEST & LOOK AT COEFFS (plot vs h): 
      np    = 100; 
      yt   = Lims(3)+(0:np)'/np*(Lims(4)-Lims(3));
      fap1  = OP_chebinterp2d(x,yt,Fcoeffs,Lims);
      %%
      for j=1:nx+1
         disp(['x = ',num2str(x(j))]);
         %%
         plot(y,Fvals(j,:).','xk'), hold on;
         %%
         plot(yt,fap1(j,:)','r'), hold off;
         pause;
      end
   end
end
