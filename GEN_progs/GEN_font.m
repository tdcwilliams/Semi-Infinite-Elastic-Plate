function GEN_font(H,fs)
%% CALL: GEN_font(H,fs);
%% H is a handle, fs is the font size, default font is 'times';

if ~exist('fs')
   fs = 16;
end
set(H,'FontName','times','FontSize',fs);
