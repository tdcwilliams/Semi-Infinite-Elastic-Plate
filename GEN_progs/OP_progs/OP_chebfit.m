function an = OP_chebfit(N,x,y);
%% CALL: an = OP_chebfit(x,y,N);
%% x is a vector of points, and the columns of y
%%  are vectors of function values
%%  corresponding to the elements of x;
%% This function creates the best least-squares fit
%%  of a Chebyshev series of degree N to the data;
%% It is equivalent to polyfit.m, but sts Chebyshev poly's are
%%  better to work with;
%% TEST: don't enter x,y and test on a trial function;

DO_TEST  = 0;
if nargin==0%% test inputs;
   DO_TEST  = 1;
   if 1
      dir0  = '~/Dropbox/MATHS/SHARED-WORK/WIM/pub_wim2d/figures/ice_charts_2/';
      fil   = [dir0,'mizwGRLFonly.mat'];
      load(fil);
      %%
      A     = mizwGRLFonly;
      %M     = length(A);
      M     = round(.6*length(A));
      time  = zeros(M,1);
      data  = zeros(M,1);
      for j=1:M
         time(j)  = A(j).date;
         data(j)  = A(j).width;
      end
      x  = time;
      y  = data;
      N  = 50;
   elseif 1
      a  = pi/4;
      b  = 3*pi/4;
      np = 100;
      x  = a+(0:np)'/np*(b-a);
      y  = sin(4*pi*x);%% works!
      %y  = x.^3-x;%% works!!
      %N  = 10;
   else%% works!!
      a  = -1;
      b  = 1;
      np = 50;
      x  = a+(0:np)'/np*(b-a);
      nt = 3;
      y  = cos(nt*acos(x));
      %N  = 10;
   end
end

if ~exist('N')
   N  = length(x)-1;
end

if (N+1)>length(x)
   disp('error in OP_chebfit.m: increase no of points or lower degree of polynomial');
   disp([N+1,length(x)]);
   disp('****');
end

dx_av = mean(x(2:end)-x(1:end-1)); 
a     = min(x)-dx_av;
b     = max(x)+dx_av;
tt    = -1+2*(x-a)/(b-a);

opt         = 0;
[TnVals,hn] = OP_interp_chebyshev(tt,{N});

if opt==1%%least squares fit;
   an = (TnVals'*TnVals)\(TnVals'*y);
else%%approx inner products;
   end_pts  = [-1;(tt(2:end)+tt(1:end-1))/2;1];
   wj       = end_pts(2:end)-end_pts(1:end-1);
   %%
   DDi   = diag(1./sqrt(hn));
   an    = DDi*(TnVals*DDi)'*(wj.*y./sqrt(1-tt.^2));
end

if DO_TEST
   plot(x,y,'-k'), hold on;
   %%
   yap      = TnVals*an;
   plot(x,yap,'r');
   %%
   %[p,s,mu] = polyfit(x,y,N);
   %xhat     = (x-mu(1))/mu(2);
   %%%
   %yap2     = p(end);
   %xx       = 1;
   %for j=1:N
   %   xx    = xx.*xhat;
   %   yap2  = yap2+p(N+1-j)*xx;
   %end
   %plot(x,yap2,'--g');
   hold off;

   [an(1)/2,mean(yap),mean(y)]
   TnVals(1:5,1)
end
