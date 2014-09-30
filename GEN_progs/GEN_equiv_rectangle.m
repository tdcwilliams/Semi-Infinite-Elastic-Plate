function [W,L,test_area,test_perimeter] = GEN_equiv_rectangle(xy,DO_PROJ)

if ~exist('DO_PROJ')
   DO_PROJ  = 0;
end

if ~exist('xy');
   if 0
      N  = 40;
      th = (1:N)'/N*2*pi;
      xy = [cos(th),sin(th)];
   else
      x  = [0 1 1 0]'*2;
      y  = [0 0 1 1]'*3;
      xy = [x y];
   end
end

if DO_PROJ==1
   %% stereographic projection from south pole;
   %% TODO: change to any centre;
   ll    = xy;
   [x,y] = GEN_stereo_proj(ll(:,1),ll(:,2));
   xy    = [x,y];
end

A  = GEN_area(xy)
P  = GEN_perimeter(xy,1)

if P^2>=16*A
   W  = (P-sqrt(P^2-16*A))/4.0;
   L  = (P+sqrt(P^2-16*A))/4.0;
   test_area      = [A,A];
   test_perimeter = [P,P];
else
   cc = [1,0,(8-A),-2*P];
   rr = roots(cc)
   W  = 0.0;
   L  = 0.0;
   for s=1:3
      a  = rr(s);
      if imag(a)==0
         W  = real(a);
         L  = W;
      end
      test_area      = [A,L*W];
      test_perimeter = [P,2*(L+W)];
   end
end
