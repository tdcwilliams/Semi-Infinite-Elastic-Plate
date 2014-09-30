% function fn_fullscreen(hd)
%
% sets the figure to full screen

function fn_fullscreen(hd)

if ~exist('hd','var')
 set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
else
 set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
end

return