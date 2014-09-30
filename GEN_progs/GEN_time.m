function t0 = GEN_time(smh)
%% CALL: GEN_time(smh)

if ~exist('smh')
   smh   = 0;
end

if smh==0%%output time taken in seconds;
   t0    = 24*60^2*rem(now,1);
elseif smh==1%%output time taken in minutes;
   t0    = 24*60*rem(now,1);
else%%output time taken in seconds;
   t0    = 24*rem(now,1);
end
