function new_record_script(file)

file=['../',file];
Directory = pwd;
if ismac
    idx = strfind(Directory,'/');
else
    idx = strfind(Directory,'\');
end

foldername = Directory(idx(end-1)+1:idx(end)-1);
if ~ismac
    [status,sheets] = xlsfinfo(file);
    nsheets = size(sheets,2);
else
    [status,cmdout] = system(['in2csv -n ', file]);
    sheets=cellstr(splitlines(convertCharsToStrings(cmdout)));
    nsheets = size(sheets,1)-1;
    sheets=sheets';
end
for i=1:nsheets
    t = datetime('today','Format','ddMMyy');
    formatOut = 'ddmmyy';
    
    sheetname = [sheets{1,i},'_',datestr(t,formatOut)];
    [~,~,RAW] = xlsread(file,sheets{1,i});
    if ~ismac
        xlswrite(file,RAW,sheetname);
    else
        xlswrite_MACOS(file,RAW,sheetname);
    end
    [RAW{:, :}]=deal(NaN);
    if ~ismac
        xlswrite(file,RAW,sheets{1,i});
    else
        xlswrite_MACOS(file,RAW,sheets{1,i});
    end
end

A={'**','ESTIMATED PARAMS', [],[],'ESTIMATED PARAMS';
    [],[],[],[],[];
    [],[],[],[],[];
    [],'EPS_GAPOILTREND_RoW',[],[],'EPS_GAPOILTREND_RoW'};
if ~ismac
    xlswrite(file,A,'All parameters & shocks');
else
    xlswrite_MACOS(file,A,'All parameters & shocks');
end
if strcmp(foldername,'GM2_commodities')
    A={ [],[],[],[],[];
        [],[],[],[],[];
        [],'1-step',[],[],[];
        'log-dens',[],[],[],'log-dens';
        [],'GYOBS_EA',[],[],'GYOBS_EA';
        [],'GC_EA',[],[],'GC_EA';
        [],'GI_EA',[],[],'GI_EA';
        [],'GG_EA',[],[],'GG_EA';
        [],'GIG_EA',[],[],'GIG_EA';
        [],'GBG_EA',[],[],'GBG_EA';
        [],'GMTOT_EA',[],[],'GMTOT_EA';
        [],'GX_EA',[],[],'GX_EA';
        [],'GL_EA',[],[],'GL_EA';
        [],'GN_EA',[],[],'GN_EA';
        [],'INOM_EA',[],[],'INOM_EA';
        [],'PHIYOBS_EA',[],[],'PHIYOBS_EA';
        [],'PHICVAT_EA',[],[],'PHICVAT_EA';
        [],'PHII_EA',[],[],'PHII_EA';
        [],'PHIG_EA',[],[],'PHIG_EA';
        [],'PHIOIL_EA',[],[],'PHIOIL_EA';
        [],'PHIM_EA',[],[],'PHIM_EA';
        [],'PHIMTOT_EA',[],[],'PHIMTOT_EA';
        [],'PHIX_EA',[],[],'PHIX_EA';
        [],'PHIW_EA',[],[],'PHIW_EA';
        [],'PHIWR_EA',[],[],'PHIWR_EA';
        [],'GYOBS_RoW',[],[],'GYOBS_RoW';
        [],'PHIYOBS_RoW',[],[],'PHIYOBS_RoW';
        [],'INOM_RoW',[],[],'INOM_RoW';
        [],[],[],[],[];
        [],'4-step',[],[],[];
        [],'GYOBS_EA',[],[],'GYOBS_EA';
        [],'GC_EA',[],[],'GC_EA';
        [],'GI_EA',[],[],'GI_EA';
        [],'GG_EA',[],[],'GG_EA';
        [],'GIG_EA',[],[],'GIG_EA';
        [],'GBG_EA',[],[],'GBG_EA';
        [],'GMTOT_EA',[],[],'GMTOT_EA';
        [],'GX_EA',[],[],'GX_EA';
        [],'GL_EA',[],[],'GL_EA';
        [],'GN_EA',[],[],'GN_EA';
        [],'INOM_EA',[],[],'INOM_EA';
        [],'PHIYOBS_EA',[],[],'PHIYOBS_EA';
        [],'PHICVAT_EA',[],[],'PHICVAT_EA';
        [],'PHII_EA',[],[],'PHII_EA';
        [],'PHIG_EA',[],[],'PHIG_EA';
        [],'PHIOIL_EA',[],[],'PHIOIL_EA';
        [],'PHIM_EA',[],[],'PHIM_EA';
        [],'PHIMTOT_EA',[],[],'PHIMTOT_EA';
        [],'PHIX_EA',[],[],'PHIX_EA';
        [],'PHIW_EA',[],[],'PHIW_EA';
        [],'PHIWR_EA',[],[],'PHIWR_EA';
        [],'GYOBS_RoW',[],[],'GYOBS_RoW';
        [],'PHIYOBS_RoW',[],[],'PHIYOBS_RoW';
        [],'INOM_RoW',[],[],'INOM_RoW';
        [],[],[],[],[];
        [],'r2',[],[],[];
        [],'GYOBS_EA',[],[],'GYOBS_EA';
        [],'GC_EA',[],[],'GC_EA';
        [],'GI_EA',[],[],'GI_EA';
        [],'GG_EA',[],[],'GG_EA';
        [],'GIG_EA',[],[],'GIG_EA';
        [],'GBG_EA',[],[],'GBG_EA';
        [],'GMTOT_EA',[],[],'GMTOT_EA';
        [],'GX_EA',[],[],'GX_EA';
        [],'GL_EA',[],[],'GL_EA';
        [],'GN_EA',[],[],'GN_EA';
        [],'INOM_EA',[],[],'INOM_EA';
        [],'PHIYOBS_EA',[],[],'PHIYOBS_EA';
        [],'PHICVAT_EA',[],[],'PHICVAT_EA';
        [],'PHII_EA',[],[],'PHII_EA';
        [],'PHIG_EA',[],[],'PHIG_EA';
        [],'PHIOIL_EA',[],[],'PHIOIL_EA';
        [],'PHIM_EA',[],[],'PHIM_EA';
        [],'PHIMTOT_EA',[],[],'PHIMTOT_EA';
        [],'PHIX_EA',[],[],'PHIX_EA';
        [],'PHIW_EA',[],[],'PHIW_EA';
        [],'PHIWR_EA',[],[],'PHIWR_EA';
        [],'GYOBS_RoW',[],[],'GYOBS_RoW';
        [],'PHIYOBS_RoW',[],[],'PHIYOBS_RoW';
        [],'INOM_RoW',[],[],'INOM_RoW';
        };
end
if strcmp(foldername,'GM3_commodities')
    
end

if strcmp(foldername,'GM2')
    
end

if strcmp(foldername,'GM3')
    
end

if strcmp(foldername,'GM3-DE')
    
end

if strcmp(foldername,'GM3-FR')
    
end

if strcmp(foldername,'GM3-IT')
    
end

if strcmp(foldername,'GM3-ES')
    
end

if strcmp(foldername,'GM3-NL')
    
end

if strcmp(foldername,'GM1-DE')
    
end

if strcmp(foldername,'GM1-FR')
    
end

if strcmp(foldername,'GM1-IT')
    
end

if strcmp(foldername,'GM1-ES')
    
end

if strcmp(foldername,'GM1-NL')
    
end
if ~ismac
    xlswrite(file,A,'RMSE''s');
else
    xlswrite_MACOS(file,A,'RMSE''s');
end

end

% fileID = fopen('All parameters and shocks.csv','w');
%
% A={'**','ESTIMATED PARAMS', [],[],'ESTIMATED PARAMS';
%     [],[],[],[],[];
%     [],[],[],[],[];
%     [],'EPS_GAPOILTREND_RoW',[],[],'EPS_GAPOILTREND_RoW'};
% [nrows, ncols] = size(A);
%
% for i = 1: nrows
%     formatSpec =[];
%     for j = 1:ncols
%         if isnumeric(A{i,j})
%             if ismissing(A{i,j})
%                 A{i,j} = [];
%             end
%             formatSpec =[formatSpec, '%9.6g,'];
%
%         elseif isstring(A{i,j}) || ischar(A{i,j})
%             formatSpec =[formatSpec, '%s,'];
%         else
%             error(message('input must be either a string or numeric'))
%         end
%     end
%
%     formatSpec =[formatSpec, '\n'];
%
%     fprintf(fileID,formatSpec,A{i,:});
% end
%
% fclose(fileID);
%
% fileID = fopen('RMSEs.csv','w');
%
% A={ [],[],[],[],[];
%     [],[],[],[],[];
%     [],'1-step',[],[],[];
%     'log-dens',[],[],[],'log-dens';
%     [],'GYOBS_EA',[],[],'GYOBS_EA';
%     [],'GC_EA',[],[],'GC_EA';
%     [],'GI_EA',[],[],'GI_EA';
%     [],'GG_EA',[],[],'GG_EA';
%     [],'GIG_EA',[],[],'GIG_EA';
%     [],'GBG_EA',[],[],'GBG_EA';
%     [],'GMTOT_EA',[],[],'GMTOT_EA';
%     [],'GX_EA',[],[],'GX_EA';
%     [],'GL_EA',[],[],'GL_EA';
%     [],'GN_EA',[],[],'GN_EA';
%     [],'INOM_EA',[],[],'INOM_EA';
%     [],'PHIYOBS_EA',[],[],'PHIYOBS_EA';
%     [],'PHICVAT_EA',[],[],'PHICVAT_EA';
%     [],'PHII_EA',[],[],'PHII_EA';
%     [],'PHIG_EA',[],[],'PHIG_EA';
%     [],'PHIOIL_EA',[],[],'PHIOIL_EA';
%     [],'PHIM_EA',[],[],'PHIM_EA';
%     [],'PHIMTOT_EA',[],[],'PHIMTOT_EA';
%     [],'PHIX_EA',[],[],'PHIX_EA';
%     [],'PHIW_EA',[],[],'PHIW_EA';
%     [],'PHIWR_EA',[],[],'PHIWR_EA';
%     [],'GYOBS_RoW',[],[],'GYOBS_RoW';
%     [],'PHYOBS_RoW',[],[],'PHIYOBS_RoW';
%     [],'INOM_RoW',[],[],'INOM_RoW';
%     [],[],[],[],[];
%     [],'4-step',[],[],[];
%     [],'GYOBS_EA',[],[],'GYOBS_EA';
%     [],'GC_EA',[],[],'GC_EA';
%     [],'GI_EA',[],[],'GI_EA';
%     [],'GG_EA',[],[],'GG_EA';
%     [],'GIG_EA',[],[],'GIG_EA';
%     [],'GBG_EA',[],[],'GBG_EA';
%     [],'GMTOT_EA',[],[],'GMTOT_EA';
%     [],'GX_EA',[],[],'GX_EA';
%     [],'GL_EA',[],[],'GL_EA';
%     [],'GN_EA',[],[],'GN_EA';
%     [],'INOM_EA',[],[],'INOM_EA';
%     [],'PHIYOBS_EA',[],[],'PHIYOBS_EA';
%     [],'PHICVAT_EA',[],[],'PHICVAT_EA';
%     [],'PHII_EA',[],[],'PHII_EA';
%     [],'PHIG_EA',[],[],'PHIG_EA';
%     [],'PHIOIL_EA',[],[],'PHIOIL_EA';
%     [],'PHIM_EA',[],[],'PHIM_EA';
%     [],'PHIMTOT_EA',[],[],'PHIMTOT_EA';
%     [],'PHIX_EA',[],[],'PHIX_EA';
%     [],'PHIW_EA',[],[],'PHIW_EA';
%     [],'PHIWR_EA',[],[],'PHIWR_EA';
%     [],'GYOBS_RoW',[],[],'GYOBS_RoW';
%     [],'PHYOBS_RoW',[],[],'PHIYOBS_RoW';
%     [],'INOM_RoW',[],[],'INOM_RoW';
%     [],[],[],[],[];
%     [],'r2',[],[],[];
%     [],'GYOBS_EA',[],[],'GYOBS_EA';
%     [],'GC_EA',[],[],'GC_EA';
%     [],'GI_EA',[],[],'GI_EA';
%     [],'GG_EA',[],[],'GG_EA';
%     [],'GIG_EA',[],[],'GIG_EA';
%     [],'GBG_EA',[],[],'GBG_EA';
%     [],'GMTOT_EA',[],[],'GMTOT_EA';
%     [],'GX_EA',[],[],'GX_EA';
%     [],'GL_EA',[],[],'GL_EA';
%     [],'GN_EA',[],[],'GN_EA';
%     [],'INOM_EA',[],[],'INOM_EA';
%     [],'PHIYOBS_EA',[],[],'PHIYOBS_EA';
%     [],'PHICVAT_EA',[],[],'PHICVAT_EA';
%     [],'PHII_EA',[],[],'PHII_EA';
%     [],'PHIG_EA',[],[],'PHIG_EA';
%     [],'PHIOIL_EA',[],[],'PHIOIL_EA';
%     [],'PHIM_EA',[],[],'PHIM_EA';
%     [],'PHIMTOT_EA',[],[],'PHIMTOT_EA';
%     [],'PHIX_EA',[],[],'PHIX_EA';
%     [],'PHIW_EA',[],[],'PHIW_EA';
%     [],'PHIWR_EA',[],[],'PHIWR_EA';
%     [],'GYOBS_RoW',[],[],'GYOBS_RoW';
%     [],'PHYOBS_RoW',[],[],'PHIYOBS_RoW';
%     [],'INOM_RoW',[],[],'INOM_RoW';
%     };
% [nrows, ncols] = size(A);
%
% for i = 1: nrows
%     formatSpec =[];
%     for j = 1:ncols
%         if isnumeric(A{i,j})
%             if ismissing(A{i,j})
%                 A{i,j} = [];
%             end
%             formatSpec =[formatSpec, '%9.6g,'];
%
%         elseif isstring(A{i,j}) || ischar(A{i,j})
%             formatSpec =[formatSpec, '%s,'];
%         else
%             error(message('input must be either a string or numeric'))
%         end
%     end
%
%     formatSpec =[formatSpec, '\n'];
%
%     fprintf(fileID,formatSpec,A{i,:});
% end
%
% fclose(fileID);