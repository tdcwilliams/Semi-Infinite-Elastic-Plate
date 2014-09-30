%% GEN_shortest_dist2line.m
%% Author: Timothy Williams
%% Date:   20130619, 13:08:32 CEST
function R1 = GEN_shortest_dist2line(z,z1,z2,varargin);

if nargin==0
   z0 = 1-1i; 
   z  = z0*[-2:2];
   z1 = 1+1i-z0;
   z2 = 0-z0;
   %R1 = [0;1;2]*abs(z0)
end

R2 = abs(z-z1);
R3 = abs(z-z2);

x1 = real(z1);
y1 = imag(z1);
x2 = real(z2);
y2 = imag(z2);
dy = y2-y1;
dx = x2-x1;

%% normal vector - eqn of line is n'*x=C
nvec  = [-dy;dx];
nvec  = nvec/sqrt(dx^2+dy^2);%%unit normal

x  = real(z);
y  = imag(z);
R1 = nvec'*[x-x1;y-y1];


if nargin==4
   %% calc distance between z and the FINITE segment between z1,z2

   %% direction along line:
   t1    = 0;
   t2    = abs(z2-z1);
   tvec  = [dx;dy];
   tvec  = tvec/sqrt(dx^2+dy^2);%%unit tangent
   %%
   dd       = [1 1i]*nvec;%%ok
   z_proj   = z-R1*dd;
   tz       = tvec'*[real(z_proj)-x1;imag(z_proj)-y1];
   %%
   if (tz>t2)|(tz<0)
      R1 = min(R2,R3);
   else
      R1 = abs(R1)
   end
else
   R1 = abs(R1);
end

if 0
   %%do some tests;
   nvec
   R1
   plot(z,'^m');
   hold on;
   plot(z_proj,'vm');
   plot([z1 z2]);
   plot([z1+R1*dd z1],'g')
   plot([z1+R1*1i*dd z1],'g')
   plot([z_proj z2],'g')
   hold off;
   daspect([1 1 1]);
end
