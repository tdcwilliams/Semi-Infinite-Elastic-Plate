function s=GEN_perimeter(xy,is_poly)

z  = xy(:,1)+1i*xy(:,2);
N  = length(z);
%%
if is_poly==0
   s  = 0*z;
   runsum=0;
   for j=2:N
     runsum=runsum+abs(z(j)-z(j-1));
     s(j)=runsum;
   end
else
   s  = abs(z(1)-z(N));
   for j=2:N
     s   = s+abs(z(j)-z(j-1));
   end
end
