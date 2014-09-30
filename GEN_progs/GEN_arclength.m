function L=GEN_arclength(rvecs)

sz=size(rvecs);
if sz(1)>sz(2)%%want points as column vectors:
  rvecs=rvecs';
  sz=fliplr(sz);
end

L=0;
for j=2:sz(2)
  dx=rvecs(1,j)-rvecs(1,j-1);
  dy=rvecs(2,j)-rvecs(2,j-1);
  L=L+abs(dx+i*dy);
end