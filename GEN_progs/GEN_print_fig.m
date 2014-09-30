function GEN_print_fig(varargin)
%% CALL: GEN_print_fig(varargin)
%% Prints (no input argument) or print-previews (1 input argument)
%% the current figure ('gcf')
saveas(gcf,'fig.eps');
if nargin==0
	!lpr -Pmasm1 fig.eps
	!rm fig.eps
else
	!gv fig.eps
end
