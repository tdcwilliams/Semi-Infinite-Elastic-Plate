function GEN_proc_fig3d(campos,cmap)

if nargin==2%% change color scheme:
  %% cmap='bone' is good
  colormap(cmap);
end

shading('interp');
set(gca,'Visible','off',...
  'CameraPosition',campos,...
    'CameraTarget',-campos);
