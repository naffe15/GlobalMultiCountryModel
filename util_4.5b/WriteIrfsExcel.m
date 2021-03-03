function []=WriteIrfsExcel(M_, oo_, xlsfilename)
% function []=WriteIrfsExcel(M_, oo_, xlsfilename)
% M_  - Dynare model M_ structure
% oo_ - Dynare output oo_ structure 
% 
if nargin<3,
    xlsfilename=[M_.fname '_Irfs.xls'];
end
if ~exist('oo_') || ~exist('M_') 
    error('quest3hlmr_results.mat')
end
% Lets go through all shocks and the irfs
irfsnames=fieldnames(oo_.irfs);
exo_names=M_.exo_names;
sheetnum=1;  % Excel sheet number 
for i=1:size(exo_names,1)
    a=strfind(irfsnames,strcat(exo_names(i,:)));
    endo_names=[];n=1;
    for j=1:size(a,1)
        if ~isempty(cell2mat(a(j)))
            k=cell2mat(a(j));
            b=strcat(irfsnames{j});
            if strcmp(b(end-length(strcat(exo_names(i,:)))+1:end),strcat(exo_names(i,:)))
                endo_names{n}=b(1:k-2);
                j=j+1;n=n+1;
            end
        end
    end
    % Write excel sheet 
    if ~isempty(endo_names)
        a=[];
        for j=1:length(endo_names)
            % sprintf('oo_.irfs.%s_%s',endo_names{j},exo_names(i,:))
            a=[a;eval(sprintf('oo_.irfs.%s_%s',endo_names{j},exo_names(i,:)))];
        end
        data=[endo_names;num2cell(a')];
        if ~ismac
            xlswrite(xlsfilename,data,exo_names(i,:)); %Write data
        else
            xlswrite_MACOS(xlsfilename,data,exo_names(i,:)); %Write data
        end
        sheetnum=sheetnum+1;
    end
end

