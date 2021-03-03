function duplicate_dynare_project(dirname)
% duplicate_dynare_project(dirname)

if nargin==0,
    disp('duplicate_dynare_project(NAME_of_NEW_DYNARE_DIR)')
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

else
    disp('new directory already exists!')
    disp('duplication failed!')
    return
end

