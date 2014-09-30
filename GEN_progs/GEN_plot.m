function GEN_plot(x,Y,c);
%% GEN_plot(x,Y,c);

for j=1:length(c)
  plot(x,Y(:,j),c{j}), hold on;
end
hold off;