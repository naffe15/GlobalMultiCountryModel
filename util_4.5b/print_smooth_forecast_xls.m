function print_smooth_forecast_xls(m_, oo_, forecasted)



fnam = fieldnames(oo_.SmoothedVariables);
MaxSize=length(fnam);

if nargin <3
    forecasted=oo_.SmoothedVariables;
    if MaxSize<250
        fname=strcat(m_.fname, '_smoothed_variables.xls');
    else
        fname1=strcat(m_.fname, '_smoothed_variables_1.xls');
        fname2=strcat(m_.fname, '_smoothed_variables_2.xls');
    end
else
    if MaxSize<250
        fname=strcat(m_.fname, '_smoothed_and_forecasted_variables.xls');
    else
        fname1=strcat(m_.fname, '_smoothed_and_forecasted_variables_1.xls');
        fname2=strcat(m_.fname, '_smoothed_and_forecasted_variables_2.xls');
    end
end



for j=1:length(fnam)
    
    raw{1,j}=fnam{j};
    
    forecastedtemp= eval(['forecasted.',fnam{j}]);
    rawvalue(:, j)=forecastedtemp;
end

if  MaxSize<250
    [success,message]=xlswrite(fname,raw,'forecasted variables');
    [success,message]=xlswrite(fname,rawvalue,'forecasted variables', 'a2');
else
    raw1=raw(1:250);
    raw2=raw(251:end);
    rawvalue1=rawvalue(:,1:250);
    rawvalue2=rawvalue(:,251:end);
    [success,message]=xlswrite(fname1,raw1,'forecasted variables');
    [success,message]=xlswrite(fname1,rawvalue1,'forecasted variables', 'a2');
    [success,message]=xlswrite(fname2,raw2,'forecasted variables');
    [success,message]=xlswrite(fname2,rawvalue2,'forecasted variables', 'a2');
end
