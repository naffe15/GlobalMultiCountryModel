function  [status] = xls_export_ann_forecasts(filename, AnnForecasts, base_name)

if nargin<3
    base_name = 'EcFin';
end
if ismac
    delete([filename '.xls'])
    delete([filename '.xlsx'])
end
    NbObs = length(AnnForecasts(1).TimeLineQA);
    NbFields = numel(fieldnames(AnnForecasts(1)));
    status = nan(length(AnnForecasts),1);
    for i=1:length(AnnForecasts)
        % Prepare the cell array to export
        cExportData = cell(NbFields, NbObs+1);
        cExportData{1,1} = 'Time';   cExportData(1,2:NbObs+1) = num2cell(AnnForecasts(i).TimeLineQA);
        if isfield(AnnForecasts(i),'EC_plot')
            cExportData{2,1} = base_name;  cExportData(2,2:NbObs+1) = num2cell(AnnForecasts(i).EC_plot);
        end
        cExportData{3,1} = 'Model';  cExportData(3,2:NbObs+1) = num2cell(AnnForecasts(i).Mean);
        cExportData{4,1} = 'HPDsup'; cExportData(4,2:NbObs+1) = num2cell(AnnForecasts(i).HPDsup);
        cExportData{5,1} = 'HPDinf'; cExportData(5,2:NbObs+1) = num2cell(AnnForecasts(i).HPDinf);
        if isfield(AnnForecasts(i),'exogassmpt')
            cExportData{6,1} = 'ExogAssmpt';  cExportData(6,2:NbObs+1) = num2cell(AnnForecasts(i).exogassmpt);
        end
        if isfield(AnnForecasts(i),'assmptrangeplot')
            xls_name = fieldnames(AnnForecasts(i).assmptrangeplot);
            for ifield=1:length(xls_name)
                cExportData{6+ifield,1} = xls_name{ifield};  cExportData(6+ifield,2:NbObs+1) = num2cell(AnnForecasts(i).assmptrangeplot.(xls_name{ifield}));
%             if isfield(AnnForecasts(i).assmptrangeplot,'High')
%                 cExportData{7,1} = 'AssmptHigh';  cExportData(7,2:NbObs+1) = num2cell(AnnForecasts(i).assmptrangeplotHigh);
%             end
%                 if
%                     isfield(AnnForecasts(i).assmptrangeplot,'Low')
%                     cExportData{8,1} = 'AssmptLow';  cExportData(8,2:NbObs+1) = num2cell(AnnForecasts(i).assmptrangeplotLow);
%                 else
%                     
%                     cExportData{7,1} = 'AssmptRange';  cExportData(7,2:NbObs+1) = num2cell(AnnForecasts(i).assmptrangeplot);
%                 end
%             end
            end
        end
        % write the sheet for a forecasted variable
        if ~ismac
            status(i)  = xlswrite(filename, cExportData, AnnForecasts(i).VarName);
        else
            writetable(array2table(cExportData), [filename '.xls'], 'Sheet', AnnForecasts(i).VarName);
%             writecell(cExportData, [filename '.xls'], 'Sheet', AnnForecasts(i).VarName );
%             status(i)  = xlswrite_MACOS(filename, cExportData, AnnForecasts(i).VarName);
        end
    end
    
end


