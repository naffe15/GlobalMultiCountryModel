function  [status] = xls_forecast_table(filename, AnnForecasts, base_name, long, sheet_name)

if nargin<3
    base_name = 'EcFin';
end
if nargin<4
   long = 0;
end
if nargin<5
   sheet_name = 'Model';
   if ismac
       delete([filename '.xls'])
       delete([filename '.xlsx'])
   end
end
diff=0;
if ~isequal(base_name , 'EcFin')
    diff=1;
end
   
xlist = {'LongName','VarName', 'TimeLineQA', 'HPDsup', 'HPDinf','assmptrangeplot'};
if long
    AnnForecasts(1).TimeLineQA =  [AnnForecasts(1).TimeLineQA(1:6); AnnForecasts(1).TimeLineQA(11);AnnForecasts(1).TimeLineQA(21);AnnForecasts(1).TimeLineQA(51);AnnForecasts(1).TimeLineQA(101)];
end
NbObs = length(AnnForecasts(1).TimeLineQA);
myFields = fieldnames(AnnForecasts(1));
myFields=myFields(~ismember(myFields,xlist));
NbFields = numel(myFields);
if long
for j=1:NbFields
    for s =1: length(AnnForecasts)
       AnnForecasts(s).(myFields{j}) =   [AnnForecasts(s).(myFields{j})(1:6); AnnForecasts(s).(myFields{j})(11); AnnForecasts(s).(myFields{j})(21); AnnForecasts(s).(myFields{j})(51); AnnForecasts(s).(myFields{j})(101)];
    end
end
end
NbVar = length(AnnForecasts);
status = nan(length(AnnForecasts),1);
mySheet=myFields;
cBaseData = {};
for k=1:NbFields
    cExportData = cell(NbVar+1, (NbObs+1));
    cExportData{1,1} = 'Time';
    cExportData(1,2:NbObs+1) = num2cell(AnnForecasts(1).TimeLineQA);
    if isequal(myFields{k},'EC_plot')
        mySheet{k}=base_name;
        
    elseif isequal(myFields{k},'Mean')
        mySheet{k}=sheet_name;
    end
    
    for i=1:length(AnnForecasts)
         if ~ismember(fieldnames(AnnForecasts(1)),'LongName')
        cExportData{i+1,1} = AnnForecasts(i).VarName;
         else
        cExportData{i+1,1} = AnnForecasts(i).LongName;
         cExportData{i+1,NbObs+2} = AnnForecasts(i).VarName;
         end
        % Prepare the cell array to export
        cExportData(i+1,2:NbObs+1) = num2cell(AnnForecasts(i).(myFields{k})-(~isequal(myFields{k},'EC_plot') && diff)*AnnForecasts(i).EC_plot);
    end
    if ~isequal(myFields{k},'EC_plot') || diff==0
        
        if ~ismac
            status(k)  = xlswrite(filename, cExportData, mySheet{k});
        else
           writetable(array2table(cExportData), [filename '.xls'], 'Sheet', mySheet{k}); 
%             writecell(cExportData,[filename '.xls'], 'Sheet', mySheet{k});
        end
    else
        cBaseData = cExportData;
        cBaseData{1,1}=base_name;
        for i=2:size(cBaseData,1)
            cBaseData{i,1}='';
        end
        
    end
end

end


