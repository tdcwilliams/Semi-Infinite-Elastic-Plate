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
if nargin==1%% test inputs;
   DO_TEST  = 1;
   if 1
      a  = pi/4;
      b  = 3*pi/4;
      np = 100;
      x  = a+(0:np)'/np*(b-a);
      y  = sin(4*pi*x);%% works!
%      y  = x.^3-x;%% works!!
   else%% works!!
      a  = -1;
      b  = 1;
      np = 50;
      x  = a+(0:np)'/np*(b-a);
      nt = 3;
      y  = cos(nt*acos(x));
   end
end

if (N+1)>length(x)
   disp('error in OP_chebfit.m: increase no of points or lower degree of polynomial');
   disp([N+1,length(x)]);
   disp('****');
end

a        = min(x);
b        = max(x);
tt       = -1+2*(x-a)/(b-a);

TnVals   = OP_interp_chebyshev(tt,{N});
an       = (TnVals'*TnVals)\(TnVals'*y);

if DO_TEST
   plot(x,y,'-k'), hold on;
   %%
   yap      = TnVals*an;
   plot(x,yap,'r');
   %%
   [p,s,mu] = polyfit(x,y,N);
   xhat     = (x-mu(1))/mu(2);
   %%
   yap2     = p(end);
   xx       = 1;
   for j=1:N
      xx = xx.*xhat;
      yap2  = yap2+p(N+1-j)*xx;
   end
   plot(x,yap2,'--g');
   hold off;
end
