function make_snapshot
% function make_snapshot
% PURPOSE Makes a snapshot from the current dynare folder
% The snapshot folder is named as pwd_snap (this should not exist) 
% and after the command the folder contains all required information to run 
% the figures (just set quest3hlmr.dyn file do_the_estimation=0)
cpwd=pwd;
snapshotdir=sprintf('%s_%s',cpwd,'snap');
duplicate_dynare_project(snapshotdir);
cd(snapshotdir);
spwd=pwd;  % Check if the folder exists
if strcmp(spwd,snapshotdir)==1
    eval(sprintf('!copy %s\\m1.mat .',cpwd));
    eval(sprintf('!copy %s\\data.mat .',cpwd));
    eval(sprintf('!copy %s\\data_calib.mat .',cpwd));
    eval('load m1.mat');
%     eval('load quest3hlmr_mode.mat');
    eval('xparam1=x(:,end);');
    eval('save quest3hlmr_mode.mat xparam1 hh -append');
else
    fprintf('Error in creating snapshot to %s',snapshotdir);
end