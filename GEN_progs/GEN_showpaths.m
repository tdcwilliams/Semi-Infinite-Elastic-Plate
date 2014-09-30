D   = path;
disp(' ');

%% AUTOMATIC PATHS ARE IN DIRECTORY defdir:
if 0%%hexagon:
   defdir   = 'local';
elseif 1%%nansen:
   defdir   = 'opt';
end
j        = [0,find(D==':')];
r        = 1;
crit     = max( D(j(r)+2:j(r)+1+length(defdir))~=defdir );
if crit==0
   disp('There are no user-defined path directories.')
else
   disp('User-defined path directories:');
   while crit
     jj    = j(r)+1:j(r+1)-1;
     disp(D(jj));
     r     = r+1;
     crit  = max( D(j(r)+2:j(r)+1+length(defdir))~=defdir );
   end
end
disp(' ');

clear D defdir j r jj crit
