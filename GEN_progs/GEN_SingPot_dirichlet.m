function [aa,ww] = GEN_SingPot_dirichlet(zz,ff,ns)
%% zz should be points in an anticlockwise direction;

if 0
   nth   = 30;
   th    = (1:30)'*2*pi/nth;
   zz    = cos(th)+2i*sin(th);
   ff    = 3*cos(th).^2-6;
   %%
end

%% get centre
%% --- this point is start point for isolines and their orthogonal lines;
nth   = length(zz);
if ~exist('ns')
   ns = nth+10;
end
%centre   = mean([zz(1),zz(round(nth/2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET SINGULAR POINTS;
dz       = zz([2:end,1])-zz;
dz_out   = -2i*dz;
ww0      = zz+dz_out;

%% at the moment we have 'nth' singularities;
%% want to interp (with respect to arclength)
%% to get 'ns' singularities;
dw    = ww0([2:end,1])-ww0;
ss    = [0;cumsum(abs(dw))];
nodes = [ww0;ww0(1)];
ss2   = linspace(0,ss(end),ns)';
ww    = GEN_lin_interp1d(ss2,ss,nodes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[W,Z]    = meshgrid(ww,zz);
Log_mat  = log(abs(W-Z));
[ii,jj]  = find(abs(Log_mat)==Inf);
if ~isempty(ii)
   Log_mat(ii,:)  = []
   Log_mat(:,jj)  = []
   if prod(size(Log_mat))==0
      aa = [];
      ww = [];
      return;
   end
end
aa       = pinv(Log_mat)*ff;

if 0
   subplot(1,2,1);
   plot(zz,'ok');
   hold on;
   plot(zz(1),'.r');
   plot(zz(2),'.g');
   %dz(1),2i*dz(1)
   %plot(ww0,'xk');
   plot(ww,'xr');
   hold off;
   %%
   subplot(1,2,2);
   ss0   = [0;cumsum( abs(dz(1:end-1)) )];
   fap   = GEN_interp_SingPot(zz,ww,aa,0);
   [ff,fap]
   plot(ss0,ff,'k');
   hold on;
   plot(ss0,fap,'--r');
   hold off;
end
