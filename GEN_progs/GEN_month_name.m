%% GEN_month_name.m
%% Author: Timothy Williams
%% Date:   20130702, 13:26:55 CEST
function MM = GEN_month_name(style)

MM = {'January','February','March','April','May','June','July','August','September','October','November','December'};
if style>0
   for j=1:12
      MM{j} = MM{j}(1:style);
   end
end
