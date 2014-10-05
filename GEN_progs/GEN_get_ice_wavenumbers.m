function k_cell=GEN_get_ice_wavenumbers(period,h,H,N)

np    = length(period);
nh    = length(h);
type  = 1;

if nh>=np
   k_cell = cell(nh,np);
   inputs = zeros(nh,3);
   for j=1:np
      inputs(:,1) = h;
      inputs(:,2) = period(j);
      inputs(:,3) = H;
      k_cell(:,j) = NDice_roots(inputs,N,type);
      disp([j np]);
   end
   k_cell   = permute(k_cell,[2 1]);
else
   k_cell   = cell(np,nh);
   inputs   = zeros(np,3);
   for j=1:nh
      inputs(:,1) = h(j);
      inputs(:,2) = period;
      inputs(:,3) = H;
      k_cell(:,j) = NDice_roots(inputs,N,type);
      disp([j nh]);
   end
end
