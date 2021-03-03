function duplicate_GM_project(dirname)
% duplicate_GM_project(dirname)

if nargin==0,
    disp('duplicate_GM_project(NAME_of_NEW_DYNARE_DIR)')
    return
end

if isempty(dir(dirname))
    mkdir(dirname)
    SUCCESS=copyfile('*.mod',dirname);
    SUCCESS=copyfile('*.dyn',dirname);
    SUCCESS=copyfile('*figures.tex',dirname);
    SUCCESS=copyfile('*report.tex',dirname);
    SUCCESS=copyfile('*.m',dirname);
    SUCCESS=copyfile('*_mode.mat',dirname);
    SUCCESS=copyfile('*_mh_mode*.mat',dirname);
    SUCCESS=copyfile('*_mean.mat',dirname);
    SUCCESS=copyfile('stodata.mat',dirname);
    SUCCESS=copyfile('pdata.mat',dirname);
    SUCCESS=copyfile('data.mat',dirname);
    SUCCESS=copyfile('dataobs.mat',dirname);
    SUCCESS=copyfile('dataoil.mat',dirname);
    SUCCESS=copyfile('data_frcst.mat',dirname);
    SUCCESS=copyfile('dataobs_frcst.mat',dirname);
    SUCCESS=copyfile('src',[dirname filesep 'src']);
    if ~ismac
        SUCCESS=copyfile('gemc\metropolis\*_mh*',[dirname filesep 'gemc' filesep 'metropolis']);
    else
        mkdir([dirname filesep 'gemc' filesep 'metropolis'])
        SUCCESS=copyfile('gemc/metropolis/*_mh*',[dirname filesep 'gemc' filesep 'metropolis'], 'f');
    end
    if SUCCESS==0,
        warning('NO METROPOLIS AVAILABLE TO DUPLICATE!')
    end

else
    disp('new directory already exists!')
    disp('duplication failed!')
    return
end
cd(dirname)

