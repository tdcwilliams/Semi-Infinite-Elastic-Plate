function nameFolds = fn_FolderNames(basedir)

if ~exist('basedir','var'); basedir='.'; end

dum_dirs=dir(basedir);
isub = [dum_dirs(:).isdir];
nameFolds = ...
 {dum_dirs(isub).name}'; nameFolds(ismember(nameFolds,{'.','..'})) = [];

return