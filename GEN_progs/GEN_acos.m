function y  = GEN_acos(x)
%%acos.m doesn't work so well close to 1;
method   = 3;

DO_TEST  = 0;
if ~exist('x')
   y0 = linspace(0,pi/2,50);
   x  = cos(y0);
   subplot(1,2,1);
   plot(x,y0);
   hold on;
   DO_TEST  = 1;
end

y     = acos(x);
dmax  = .025;
%%
jj    = find(abs(1-x)<=dmax & abs(1-x)>0);
Y     = y(jj);
%%
a1 = .5;
a2 = -a1/4/3;
a3 = -a2/6/5;
%%
switch method

case 1%%near x=1, approx cos as a polynomial in x^2;
      %%solve with roots.m;
      %%solution is the closest one to the matlab solution;
   for k=1:length(jj)
      d        = 1-x(jj(k));
      Yk       = Y(k);
      if 1%%cubic (<2.5e-10 accuracy)
         cc = [a3 a2 a1 -d];
      else%%quadratic (<4e-7 accuracy)
         cc = [a2 a1 -d];
      end
      rc       = roots(cc);
      ysq      = rc(find(abs(rc-Yk.^2)==min(abs(rc-Yk.^2))));
      %[k ysq]
      y(jj(k)) = sqrt(ysq);
   end
case 2%%near x=1, approx acos as a quadratic in x^2 (<4e-7 accuracy);
      %%analytical solution;
      %%solution is the closest one to the matlab solution;
   dd    = 1-x(jj);
   ysq   = ( -a1+sqrt(a1^2+4*a2*dd) )/2/a2;
   if 0%%don't seem to need this solution;
      ysq2  = ( -a1-sqrt(a1^2+4*a2*dd) )/2/a2;
      for k=1:length(jj);
         if abs(ysq2(k)-Y(k)^2)<abs(ysq(k)-Y(k)^2)
            ysq(k)   = ysq2(k);
         end
         %[k ysq(k)]
      end
   end
   y(jj) = sqrt(ysq);
case 3%%use Muller's method (<8e-16 accuracy);
   for k=1:length(jj)
      %[k x(jj(k))]
      rr = GEN_findroots_muller( @fn_zero,Y(k),x(jj(k)) );
      %pause
      y(jj(k)) = rr;
   end
   %y(jj) = GEN_findroots_muller( @fn_zero,Y,x(jj) );
end

if DO_TEST
   plot(x,y,'r');
   hold off;
   %%
   subplot(1,2,2);
   plot(x,abs(y-y0));

   rr = [GEN_findroots_muller( @fn_zero,0,1 ) 0]
end




function z  = fn_zero(y,x)
%% near 1, cos(x) is a quadratic, so use muller's method; 
z  = cos(y)-x;
