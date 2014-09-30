function nameFolds = fn_FileNames(basedir)

if ~exist('basedir','var'); basedir='.'; end

dum_dirs=dir(basedir);
isub = ~[dum_dirs(:).isdir];
nameFolds = ...
 {dum_dirs(isub).name}';

return