function LatexPwd()

path0=pwd;
i=strfind(path0,filesep);
path0=path0(i(end)+1:end);
% path0=strrep([pwd,'/'],'\','/');
fid=fopen('pwd.tex', 'w+');
fprintf(fid,['\\usepackage[maindir=/',path0,'/]{currfile}\n']);
fclose(fid);