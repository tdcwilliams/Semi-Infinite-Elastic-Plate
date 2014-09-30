function y = GEN_SingPot_interp(zz,ww,aa,ndiff);
%% CALL: y = GEN_interp_SingPot(zz,ww,aa,ndiff);

[W,Z]    = meshgrid(ww,zz);
if ndiff==0
   Log_mat  = log(abs(W-Z));
   y        = Log_mat*(aa);
else
   R        = abs(Z-W);
   dLog_mat = (Z-W)./R.^2;
   y        = dLog_mat*(aa);%%this gives F_x+1i*F_y;
end
