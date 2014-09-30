function GEN_move_label(x_or_y,amount)

if x_or_y=='x'%%move xlabel up or down
  X=get( gca, 'xlabel');
  v=get(X,'position');
  if length(amount)==1
    r=2;
    v(r)=v(r)+amount;
  else
    v=v+amount;
  end
  set(X,'position',v);
else%%move ylabel left or right
  X=get( gca, 'ylabel');
  v=get(X,'position');
  if length(amount)==1
    r=1;
    v(r)=v(r)+amount;
  else
    v=v+amount;
  end
  set(X,'position',v);
end
