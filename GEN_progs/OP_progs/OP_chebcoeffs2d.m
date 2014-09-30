function [F_coeffs,tt,ww]=OP_chebcoeffs2d(tt,ww,Fmat0);


DO_TEST  = 0;
if nargin==0
   DO_TEST  = 1;
   nt       = 10;
   [t,w]    = OP_numint_chebyshev(nt+1);
   tt       = {t,t};
   ww       = {w,w};
   Fmat0    = sin(2*pi*t)*cos(pi*t');
end


N(1)     = size(Fmat0,1)-1;
N(2)     = size(Fmat0,2)-1;

for j=1:2
   if isempty(tt{j})
      [t,w] = OP_numint_chebyshev(N(j)+1);
      tt{j} = t;
      ww{j} = w;
   else
      t     = tt{j};
      w     = ww{j};
   end
   ipT{j}   = OP_inprod_chebyshev(t,w,N(j));
end

F_coeffs = ipT{1}*Fmat0*ipT{2}';

if DO_TEST
   np = 100;
   t  = -1+2*(0:np)'/np;
   for j=1:nt+1
      plot(tt{1},Fmat0(:,j),'.k'), hold on;
      %%
      fap   = OP_chebinterp2d(t,tt{2}(j),F_coeffs,[-1,1,-1,1]);
      plot(t,fap), hold off;
      pause;
   end
end
