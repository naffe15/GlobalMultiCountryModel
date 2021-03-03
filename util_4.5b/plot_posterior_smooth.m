function plot_posterior_smooth(varargin)
global M_ oo_

dname = 'metropolis';
% dname = 'SLICE';

varlist = char(varargin);
DirectoryName = CheckPath(dname,M_.dname);

load([DirectoryName '/' M_.fname '_data.mat'],'stock_gend','stock_data');

temp_smooth_file_list = dir([DirectoryName,filesep,'*_smooth*.mat']);
jfile=0;
for j=1:length(temp_smooth_file_list),
    if isempty(strfind(temp_smooth_file_list(j).name,'smoothed')),
        jfile=jfile+1;
        ffil(jfile)=temp_smooth_file_list(j);
    end
end
% ffil = dir([DirectoryName '/' M_.fname '_smooth*.mat']);
fileParam = dir([DirectoryName '/' M_.fname '_param*.mat']);
fileShock = dir([DirectoryName '/' M_.fname '_inno*.mat']);

x=[];
yss = [];
logpo = [];
for j=1:length(fileParam),
    try,
    load([DirectoryName '/' M_.fname '_param' int2str(j) '.mat'],'stock','stock_logpo','stock_ys');
    x = [x; stock];
    yss = [yss; stock_ys];
    logpo = [logpo; stock_logpo];
    catch
    end
end

B = size(x,1);

pm3(M_.endo_nbr,stock_gend,length(ffil),B,'Smoothed variables',...
	'',varlist,M_.endo_names_tex,M_.endo_names,...
    varlist,'SmoothedVariables',CheckPath('metropolis',M_.dname),'_smooth'); 

