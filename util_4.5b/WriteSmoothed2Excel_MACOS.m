function []=WriteSmoothed2Excel_MACOS(M_, oo_, xlsfilename)


if nargin<3,
    xlsfilename=[M_.fname '_Smoothed.xls'];
end


load data

InitialTime=T(1);
clear data


if ~exist('oo_') || ~exist('M_') 
    error('quest3hlmr_results.mat')
end



VarNames=fieldnames(oo_.SmoothedVariables);
SmoothedVariables=[struct2cell(oo_.SmoothedVariables)];
my_field_names = [fieldnames(oo_.SmoothedVariables)];
isvar=zeros(length(SmoothedVariables),1);
for jf = 1:length(SmoothedVariables),
    isvar(jf)=~(isstruct(SmoothedVariables{jf}));
end
VarNames=VarNames(logical(isvar));

TempVarValues=getfield(oo_.SmoothedVariables,VarNames{1,:});
VarLength=size(TempVarValues,1);

EndTimeVar=InitialTime+fix(VarLength/4)+(VarLength-4*fix(VarLength/4))/4;
TimeSpanVar=InitialTime:0.25:EndTimeVar;
TimeSpanVar=TimeSpanVar(1:end-1)';

NbrVars=size(VarNames,1);
VarNamesIndx=ones(NbrVars,1);
for jj=1:NbrVars
    
    if findstr('AUX_ENDO_LAG_',VarNames{jj,:})
        VarNamesIndx(jj,1)=0;
    end
end

VarValues={};
VarValues(1,1)={'Var Name'};
VarValues(2:length(TimeSpanVar)+1,1)=num2cell(TimeSpanVar)';
VarNames1=VarNames(find(VarNamesIndx==1));
for jj=1:size(VarNames1,1)
VarValues(1,jj+1)=VarNames1(jj,:);
VarValues(2:VarLength+1,jj+1)=num2cell(getfield(oo_.SmoothedVariables,VarNames1{jj,:}));

end



ShocksNames=fieldnames(oo_.SmoothedShocks);
SmoothedShocks=[struct2cell(oo_.SmoothedShocks)];
my_field_names = [fieldnames(oo_.SmoothedShocks)];
isvar=zeros(length(SmoothedShocks),1);
for jf = 1:length(SmoothedShocks),
    isvar(jf)=~(isstruct(SmoothedShocks{jf}));
end
ShocksNames=ShocksNames(logical(isvar));
TempShocksValues=getfield(oo_.SmoothedShocks,ShocksNames{1,:});
ShocksLength=size(TempShocksValues,1);
EndTimeShocks=InitialTime+fix(ShocksLength/4)+(ShocksLength-4*fix(ShocksLength/4))/4;
TimeSpanShocks=InitialTime:0.25:EndTimeShocks;
TimeSpanShocks=TimeSpanShocks(1:end-1)';



NbrShocks=size(ShocksNames,1);

ShockValues={};
ShockValues(1,1)={'Shock Name'};
ShockValues(2:length(TimeSpanShocks)+1,1)=num2cell(TimeSpanShocks)';

for jj=1:size(ShocksNames,1)
ShockValues(1,jj+1)=ShocksNames(jj,:);
ShockValues(2:VarLength+1,jj+1)=num2cell(getfield(oo_.SmoothedShocks,ShocksNames{jj,:}));
end

fileID = fopen('Variables.csv','w');
formatSpec1 = [ '%s, ', repmat('%6.2f, ',1,size(VarValues,1)-1),'\n'];
formatSpec = [ '%s, ', repmat('%9.6g, ',1,size(VarValues,1)-1), '\n'];

nrows = size(VarValues,2);
fprintf(fileID,formatSpec1,VarValues{:,1});
for ii = 2: nrows
    fprintf(fileID,formatSpec,VarValues{:,ii});
end

fclose(fileID);


fileID = fopen('Shocks.csv','w');
formatSpec1 = [ '%s, ', repmat('%6.2f, ',1,size(ShockValues,1)-1),'\n'];
formatSpec = [ '%s, ', repmat('%9.6g, ',1,size(ShockValues,1)-1), '\n'];

nrows = size(ShockValues,2);
fprintf(fileID,formatSpec1,ShockValues{:,1});
for ii = 2: nrows
    fprintf(fileID,formatSpec,ShockValues{:,ii});
end

fclose(fileID);
% tempTable= cell2table(VarValues(:,2:end)','VariableNames',(cellstr(xtemp))');
% writetable(tempTable,'Vars.csv');



%xlswrite(xlsfilename,VarValues','Variables')
%xlswrite(xlsfilename,ShockValues','Shocks')
