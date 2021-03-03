
load quest3hlmr_results

date_now = clock;
if date_now(5) < 10
    date_now = strcat(num2str(date_now(1)),'_',num2str(date_now(2)),'_', ...
        num2str(date_now(3)), '_at_',num2str(date_now(4)),'_', '0',num2str(date_now(5)));
    
else
    date_now = strcat(num2str(date_now(1)),'_',num2str(date_now(2)),'_', ...
        num2str(date_now(3)), '_at_',num2str(date_now(4)),'_', num2str(date_now(5)));
    
end


folders=pwd;
% % % % % latexfolder=strcat('Latex_', date_now);
% % % % % mkdir(latexfolder);
% % % % % 
% % % % % cd(latexfolder);
% % % % % latexdir=pwd;
% % % % % cd ..

pos=strfind(folders,'\');
dirname=folders(pos(end)+1:end);

folders1=strcat(folders,'\','*.tex');



srtlist = fuf (folders1,1,'detail');
% fnamereport=strcat(M_.fname,'_Report_','.tex');
fnamereport=strcat(dirname,date_now,'_Report','.tex');
% fullreportname=strcat(latexdir,'\',fnamereport);
fidtex=fopen(fnamereport,'w+');

fprintf(fidtex,'%s\n','\documentclass[10pt,a4paper]{article}');
fprintf(fidtex,'%s\n',' \usepackage{geometry}');
fprintf(fidtex,'%s\n','\usepackage{fullpage}');
fprintf(fidtex,'%s\n','\usepackage{psfrag,graphicx, graphics}');
fprintf(fidtex,'%s\n','\usepackage{hyperref}');
fprintf(fidtex,'%s\n','\usepackage{float}');
fprintf(fidtex,'%s\n','\usepackage{perpage}');
fprintf(fidtex,'%s\n','\MakeSorted{figure}');
fprintf(fidtex,'%s\n','\MakeSorted{table}');
fprintf(fidtex,'%s\n','\usepackage{setspace}');

fprintf(fidtex,'%s\n','\begin{document}');

fprintf(fidtex,'%s\n','\doublespacing');
fprintf(fidtex,'%s\n','\title{QUEST}');
fprintf(fidtex,'%s\n','\author{ECFIN-JRC-QUEST Team\\DG ECFIN, Joint Research Centre \\ European Commission}');

fprintf(fidtex,'%s\n','\maketitle');

fprintf(fidtex,'%s\n','\footnotesize');

fprintf(fidtex,'%s\n', '\href{quest3hlmr_dynamic.pdf }{Open the model definition file}');

% Use \titlerunning{Short Title} for an abbreviated version of
% your contribution title if the original one is too long




% copyfile('myFun.m','d:/work/Projects/')

for j=1:length(srtlist)
    strlistj=srtlist{j};
    matchrep= findstr(strlistj,'Report.tex');
    matchdyn=findstr(strlistj, 'dynamic');
    matchstat=findstr(strlistj ,'static');
    if  ~isempty(matchrep) || ~isempty(matchdyn)  || ~isempty(matchstat)
        fprintf(fidtex,'%s\n','');
    else
        [a,b,c]= fileparts(strlistj);
        includedfile0=strcat(a,'\',b) ;
        includedfile0=regexprep(includedfile0, '\','/');
        includedfile=strcat('\include{',includedfile0,'}');
        fprintf(fidtex,'%s\n','');
        fprintf(fidtex,'%s\n',includedfile);
        fprintf(fidtex,'%s\n','');
    end
end


fprintf(fidtex,'%s\n','\end{document}');

fclose(fidtex);



% open(fnamereport)
   [a,b,c]= fileparts(fnamereport);
   
[dump1, dump2]=system(['latex ',b]);
[dump1, dump2]=system(['dvips ',b]);
[dump1, dump2]=system(['ps2pdf ',strcat(b, '.ps')]);
