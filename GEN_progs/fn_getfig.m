% function fig=fn_getfig
%
% DESCRIPTION: find the next figure handle available
%
% L Bennetts Sept 2013 / Adelaide

function fig=fn_getfig(num)

if ~exist('num','var'); num=1; end

figHandles = findobj('Type','figure');

if isempty(figHandles)
 fig=1:num;
else
 dum = 1:max(figHandles)+num;
 dum(figHandles)=[];
 fig = dum(1:num);
end

return