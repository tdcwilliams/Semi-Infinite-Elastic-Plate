function GEN_colormap(IF_BW,vartype)
%% CALL: GEN_colormap(IF_BW,vartype)
%% IF_BW==1 for black and white,
%% IF_BW==0 for colour;
%% vartype==1 is good for Hs
%% vartype==2 is good for Dmax

if ~exist('vartype')
   vartype  = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IF_BW%% Black & white plots;
   if vartype==1%% Hs
      %colormap bluetone%M
      %colormap seacol2%m
      %colormap seacol3%M
      %colormap seacol_inv%m
      %%
      colormap bone%M
      %colormap gray%M
      %colormap gray%M
      %colormap copper
      %colormap jet%=ncview
      return
   end
   if vartype==2%% Dmax
      %colormap bluetone%M
      %colormap gray%M
      colormap copper
      %colormap gray%M
      %colormap summer
      return;
   end
   %colormap matsmap
   colormap bluetone%M
   %colormap seacol2
   %colormap seacol3%m
   %colormap earthc
   %colormap gaudy
   %colormap fc100
   %colormap helix
   %colormap helix2
   %colormap seacol
   %colormap seacol_inv%m
   %colormap ssec
   %colormap sstpal
   %%
   %colormap hsv
   %colormap hot
   %colormap gray%M
   %colormap copper
   %colormap vga
   %colormap jet%=ncview
   %colormap prism
   %colormap cool
   %colormap autumn
   %colormap spring
   %colormap winter
   %colormap summer
   %%colormap pink
   %%colormap white
   %%colormap flag
   %%colormap lines
   %%colormap colorcube
   %%colormap bone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else%% COLOR PLOTS
   if vartype==1%% Hs
      %colormap bluetone%M
      colormap seacol3%M
      %colormap helix
      %colormap helix2
      %colormap bone%m
      %colormap copper
      %colormap jet%=ncview
      %colormap cool
      %colormap summer
      return;
   end
   if vartype==2%% Dmax
      colormap bluetone%M
      %colormap seacol3%m
      %colormap copper
      %colormap jet%=ncview
      %colormap summer
      return;
   end
   colormap bluetone%M
end
