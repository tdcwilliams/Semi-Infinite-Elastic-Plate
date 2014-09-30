function GEN_progrep(nN,t0,smh)
%% CALL: GEN_progrep(nN,t0,smh)
%% t0,smh optional;

disp(' ');
disp([num2str(nN(1)),...
      ' runs completed (out of ',...
      num2str(nN(2)),').']);

if  nargin==1
   return;
end
if ~exist('smh')
   smh   = 0;
end

str1  = 'Time taken: ';
if smh==0%%output time taken in seconds;
   t1    = 24*60^2*rem(now,1);
   str2  = 's.';
elseif smh==1%%output time taken in minutes;
   t1    = 24*60*rem(now,1);
   str2  = ' mins.';
else%%output time taken in seconds;
   t1    = 24*rem(now,1);
   str2  = 'h.';
end

dt = round((t1-t0)*10)/10;
disp([str1,num2str(dt),str2]);
