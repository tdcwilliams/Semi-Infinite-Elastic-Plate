function GEN_man_cbar(posvec,y_rng)

subplot('position',posvec);
y     = linspace(y_rng(1),y_rng(2),20);
x     = linspace(0,1,20);
[Y,X] = meshgrid(y,x);
Z     = Y;
pcolor(X,Y,Z);
set(gca,'fontname','times','fontsize',16);
set(gca,'xticklabels',[],'yaxislocation','right');
shading interp
box on;
