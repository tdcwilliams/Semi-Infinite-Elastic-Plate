function txt = GEN_text(ax,ay,str)

xl = get(gca,'xlim');
yl = get(gca,'ylim');

x  = xl(1)+ax*(xl(2)-xl(1));
y  = yl(1)+ay*(yl(2)-yl(1));

txt   = text(x,y,str);
