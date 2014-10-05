function [an,y_approx,int_prams]  = OP_arbpts_legendre(time,data,N)

DO_TEST  = 0;
if nargin==0
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
      N  = 4;
   else
      time  = linspace(-1,1,60)';
      data  = 1./(1+.5*sin(3*time));
      N     = 11;
   end
   DO_TEST  = 1;
   %N        = 10;
end

M  = length(time);
if ~exist('N')
   N  = M-1;
end

dt    = time(2:end)-time(1:end-1);
dt_av = mean(dt);

%%normalised time interval;
a           = time(1)-dt_av/2;
b           = time(end)+dt_av/2;
T           = b-a;
int_prams   = [a b];
tj          = -1+2/T*(time-a);

IS_ORTHONORMAL = 1;%% normalise so that <P_m,P_n>=\delta_{mn}

Pm = OP_interp_legendre(tj,{N},IS_ORTHONORMAL);%%M X N+1
if 0%%least squares;
   f  = Pm'*data;%%N+1 x 1
   K  = Pm'*Pm;%%N+1 x N+1
   an = K\f;%% [an;bn]
else%%approx inner product
   end_pts  = [-1;(tj(2:end)+tj(1:end-1))/2;1];
   wj       = end_pts(2:end)-end_pts(1:end-1);
   an       = Pm'*(wj.*data);
end

if nargout>1
   y_approx = Pm*an;
end

if DO_TEST==1
   an
   %Pm(1:5,1)
   if 0%%use all coefficients
      N1 = 0;
      N2 = N;
   elseif 0%%detrend;
      N1 = 0;%%N1=0 uses all coeff's, N1=1 ignores the a0 term, N1=2 ignores a0,a1 and b1,...
      N2 = N;%%upper truncation
   else%%look at low-freq signal;
      N1 = 0;
      N2 = N;
   end
   jt1      = 1+(N1:N2)';
   jt2      = (max(N1,1):N2)';
   approx   = Pm(:,jt1)*an(jt1);
   [an(1)/2,mean(data),mean(approx)]
   mean(time)
   plot(time,data,'-k',...
        time,approx,'--r',...
        time,0*time,'-b');
end

