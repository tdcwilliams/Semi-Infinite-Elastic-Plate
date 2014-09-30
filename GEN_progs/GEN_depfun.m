function list2 = GEN_depfun(mfile)
%% CALL: GEN_depfun(mfile)
%% list non-standard file dependencies of 'mfile' and put them in
%% a directory

mydir    = '/home/nersc/timill/matlab';
n        = length(mydir);
%%
%outfile  = 'dep.sh';%%script to be created
outdir   = input('Name of folder to copy files to (enter a string): ');
eval(['!mkdir ' outdir]);

%% get dependent functions and eliminate built-in functions;
list     = depfun(mfile);
N        = length(list);
ll = 0;
for j=1:N
   ff = list{j};
   if strcmp(mydir,ff(1:n))
      %disp(ff);
      list2{j} = ff;
      ll = max(ll,length(ff));

      eval(['!cp ' ff  ' ' outdir]);
   end
end

disp(['Files should be in ' outdir ' directory,']);
disp('Files:');
eval(['!ls -lh ' outdir]);

%zz = blanks(ll+3);
%%
%disp(['This function creates a script (' outfile ')']);
%disp( 'that copies the dpendent files of');
%disp(mfile);
%disp('to the folder of your choice.');
%disp(' ');

%%% create script:
%fid      = fopen(outfile, 'w');
%%%
%s  = ['# Script to collect dependent files in ', outdir,...
%      ' and ',outdir,'.tar;'];
%fprintf(fid,[s '\n']);%%'\n' is newline
%%%
%s  = ['# Edit now, if required, before it is executed.'];
%fprintf(fid,[s '\n']);%%'\n' is newline
%fprintf(fid,['\n']);
%%%
%for j=1:length(list2)
%   s0                = zz;
%   ff                = list2{j};
%   s0(1:length(ff))  = ff;
%   s  = ['cp   ',s0,outdir];
%   fprintf(fid,[s '\n']);
%end
%
%fprintf(fid,['\n']);
%s  = ['tar -cvf ',outdir,'.tar ',outdir];
%fprintf(fid,[s]);
%fclose(fid);
%
%eval(['!vim ',outfile]);
%%%
%eval(['!chmod 750 ',outfile]);
%disp(' ');
%disp(['Files should be in ' outdir ' directory,']);
%disp('and in a tar file made with:');
%disp(s);
%disp(' ');
%%%
%disp('Contents of tar file:');
%disp(' ');
%eval(['!./',outfile]);
