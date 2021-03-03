function make_GM_slice_snapshot
% function make_GM_slice_snapshot
% PURPOSE Makes a snapshot from the current dynare folder
% The snapshot folder is named as pwd_snap (this should not exist) 
% and after the command the folder contains all required information to run 
% the figures (just set quest3hlmr.dyn file do_the_estimation=0)
cpwd=pwd;
snapshotdir=sprintf('%s_%s',cpwd,'snap');
duplicate_dynare_project(snapshotdir);
SUCCESS=copyfile('data.mat',snapshotdir);
SUCCESS=copyfile('dataobs.mat',snapshotdir);
SUCCESS=copyfile('dataoil.mat',snapshotdir);
SUCCESS=copyfile('stodata.mat',snapshotdir);
SUCCESS=copyfile('gemc',[snapshotdir filesep 'gemc']);
SUCCESS=copyfile('src',[snapshotdir filesep 'src']);
cd(snapshotdir);
spwd=pwd;  % Check if the folder exists
if strcmp(spwd,snapshotdir)==1
    
    a=dir('gemc/metropolis/gemc_mh_tmp*.*');
    b=dir('gemc/metropolis/gemc_mh1_blck*.*'); % some chain may be already complete ...
    a=[a; b];
    m=-inf;
    for ja=1:length(a)
        load(['gemc/metropolis/' a(ja).name])
        [m1, im]=max(logpo2(find(logpo2)));
        if m1>m
            xparam1=x2(im,:)';
            % update fval at mode
            m=m1;
        end
    end        
    if ~ismac
    load( [ cpwd '\gemc\prior\definition.mat']);
    else
    load( [ cpwd '/gemc/prior/definition.mat']);    
    end
    parameter_names = bayestopt_.name;
    hh=eye(length(xparam1));
    save('gemc_mh_mode','xparam1','hh','parameter_names','-append')
else
    fprintf('Error in creating snapshot to %s',snapshotdir);
end