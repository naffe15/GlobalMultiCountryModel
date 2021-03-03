function [output_file] = truncate_forecast_data(input_file,last_data,sample_size,assumptions)

% This function takes a file with confidential forecasts and replaces out
% of sample data with n.a.
% It uses as inputs the name of the file, the last date of the sample and
% the overall sample size including the forecast horizon
% If assumptions is used, it keeps the lines corresponding to the
% assumptions, that is INOM_EA, E_EA, POIL_RoW, YOBS_RoW, YOBS_RoW+US,
% E_RoW, E_RoW+US
% Example: Duplicate_insample_vintage('GM2016Q4.xlsx','2016Q4',sample_size,assumptions)

%% Check the inputs of the function

%% A) Check the format of the input_file
input_file_ext = input_file(end-4:end);
binary = strcmp(input_file_ext,'.xlsx');
if binary==0
   error('Error. Input_file must be an excel file with extension .xlsx')
end

%% B) Check the format of the last_data argument
last_data_length = length(last_data);
last_data_format = last_data(5:5);
last_data_bin = strcmp(last_data_format,'Q');
last_data_quarter = last_data(6:6);
if last_data_length ~= 6 || last_data_bin ~= 1 || str2num(last_data_quarter) < 1 || str2num(last_data_quarter) > 4
   error('Error. Last_data must be a string of the type yyyyQx where yyyy is the year and x is the quarter')
end

%% C) Check the length of the sample_size
sample_size_mod = mod(sample_size,4);
if sample_size_mod ~= 0
   error('Error. Sample_size must be a multiple of 4')
end

%% D) Handle assumptions
if nargin<4,
    assumptions=0;
end

if assumptions
    % INOM_EA
    assmpt_info(1).co = 'EA19';
    assmpt_info(1).name = 'INOM';
    assmpt_info(1).blck = 'Q';
    assmpt_info(1).row = '2';
    
    % E_EA
    assmpt_info(2).co = 'EA19';
    assmpt_info(2).name = 'E';
    assmpt_info(2).blck = 'Q';
    assmpt_info(2).row = '66';
    
    % POIL_US
    assmpt_info(3).co = 'US';
    assmpt_info(3).name = 'POIL';
    assmpt_info(3).blck = 'Q';
    assmpt_info(3).row = '48';
    
    % YOBS_RoW
    assmpt_info(4).co = 'RoW';
    assmpt_info(4).name = 'YOBS';
    assmpt_info(4).blck = 'YR';
    assmpt_info(4).row = '1';

    % YOBS_RoW+US
    assmpt_info(5).co = 'RoW+US';
    assmpt_info(5).name = 'YOBS';
    assmpt_info(5).blck = 'YR';
    assmpt_info(5).row = '1';

    % E_RoW
    assmpt_info(6).co = 'RoW';
    assmpt_info(6).name = 'E';
    assmpt_info(6).blck = 'YR';
    assmpt_info(6).row = '4';

    % E_RoW+US
    assmpt_info(7).co = 'RoW+US';
    assmpt_info(7).name = 'E';
    assmpt_info(7).blck = 'YR';
    assmpt_info(7).row = '4';
end

input_file_a = input_file(1:end-5);
output_file = strcat(input_file_a,'_insample.xlsx');
SUCCESS = copyfile(input_file,output_file);


%% First annual data block
xlRangeA='K2:GA43';
[numA,txtA,rawA] = xlsread(input_file,'EA19',xlRangeA);
last_year = 1995 + sample_size/4 - 1;
year = str2num(last_data(1:end-2));
rawA = rawA(:,1:(last_year-1995+1));
maxyear = 1995 + size(rawA,2);
quarter = str2num(last_data(end:end));


if (year - 2010)<=26
    if quarter == 4
      column = strcat('A',char(65 + (year - 2010)));
      fore_years = maxyear - year - 1;
    else
      column = strcat('A',char(65 + year - 2010 - 1));
      fore_years = maxyear - year;
    end
else
    if quarter == 4
      column = strcat('B',char(65 + year - 2036));
      fore_years = maxyear - year - 1;
    else
      column = strcat('B',char(65 + year - 2036 - 1));
      fore_years = maxyear - year;
    end
end
    
xlRange = strcat(column,'2');
clearvars x;

for i=1:size(rawA,1)
    for j=1:fore_years
        x{i,j} = 'n.a.';
    end
end

countries = {'EA19','DE','FR','IT','ES','US'};

for c=1:length(countries)
    co=countries{c};
    xlswrite(output_file,x,co,xlRange);
end

%% Quarterly data block
xlRangeQ='K52:GA121';
[numQ,txtQ,rawQ] = xlsread(input_file,'EA19',xlRangeQ);
rawQ = rawQ(:,1:sample_size);
year = str2num(last_data(1:end-2));
quarter = str2num(last_data(end:end));
fore_quarters = sample_size - (year-1995)*4 - quarter;

if ((year - 2012)*4 + quarter) <= 26
   column = strcat('C',char(65 + (year - 2012)*4 + quarter));
else
   column = strcat('D',char(65 + (year - 2018.5)*4 + quarter));
end

xlRange = strcat(column,'52');
clearvars x;

for i=1:size(rawQ,1)
    for j=1:fore_quarters
    x{i,j} = 'n.a.';
    end
end

countries = {'EA19','DE','FR','IT','ES','US'};

for c=1:length(countries)
    co=countries{c};
    xlswrite(output_file,x,co,xlRange);
end

%% Second annual data block
xlRangeA='K152:GA221';
[numA,txtA,rawA] = xlsread(input_file,'EA19',xlRangeA);
last_year = 1995 + sample_size/4 - 1;
year = str2num(last_data(1:end-2));
rawA = rawA(:,1:(last_year-1995+1));
maxyear = 1995 + size(rawA,2);
quarter = str2num(last_data(end:end));

if (year - 2010)<=26
    if quarter == 4
      column = strcat('A',char(65 + (year - 2010)));
      fore_years = maxyear - year - 1;
    else
      column = strcat('A',char(65 + year - 2010 - 1));
      fore_years = maxyear - year;
    end
else
    if quarter == 4
      column = strcat('B',char(65 + year - 2036));
      fore_years = maxyear - year - 1;
    else
      column = strcat('B',char(65 + year - 2036 - 1));
      fore_years = maxyear - year;
    end
end
    
xlRange = strcat(column,'152');
clearvars x;

for i=1:size(rawA,1)
    for j=1:fore_years
        x{i,j} = 'n.a.';
    end
end

countries = {'EA19','DE','FR','IT','ES','US'};

for c=1:length(countries)
    co=countries{c};
    xlswrite(output_file,x,co,xlRange);
end


%% RoW e RoW+US
countries = {'RoW','RoW+US'};
for c=1:length(countries)
    co=countries{c};

xlRangeRoW='A3:M35';
[numRoW,txtRoW,rawRoW] = xlsread(input_file,co,xlRangeRoW);
year = str2num(last_data(1:end-2));
%maxyear = max(numRoW(:,1));
maxyear = 1995 + size(numRoW,1) - 1;
quarter = str2num(last_data(end:end));

if quarter == 4
   fore_years = maxyear - year;
else
   fore_years = maxyear - year + 1;
end

clearvars x;

for i=1:fore_years
    for j=1:size(rawRoW,2)-1
        x{i,j} = 'n.a.';
    end
end

if quarter == 4
   row = num2str(year - 1995 + 4);
else
   row = num2str(year - 1995 + 4 - 1);
end

xlRange = strcat('B',row);
xlswrite(output_file,x,co,xlRange);
clearvars x;

for i=1:fore_years
    for j=1:size(rawRoW,2)-5
        x{i,j} = 'n.a.';
    end
end

if quarter == 4
   row = num2str(year - 1995 + 44 + 1);
else
   row = num2str(year - 1995 + 44 + 1 - 1);
end

xlRange = strcat('B',row);
xlswrite(output_file,x,co,xlRange);
clearvars x;
end

%% Writes the assumptions in the Excel file
if assumptions

for m=1:size(assmpt_info,2)
    if assmpt_info(m).blck == 'Q';
        xlRangeQ='K52:GA119';
        [numQ,txtQ,rawQ] = xlsread(input_file,assmpt_info(m).co,xlRangeQ);
        rawQ = rawQ(:,1:sample_size);
        var = rawQ(str2num(assmpt_info(m).row),:);
        line= num2str((51+str2num(assmpt_info(m).row)));
        xlRange = strcat('K',line);
        xlswrite(output_file,var,assmpt_info(m).co,xlRange);
    elseif assmpt_info(m).blck == 'YR';
        xlRangeA='B3:M33';
        [numA,txtA,rawA] = xlsread(input_file,assmpt_info(m).co,xlRangeA);
        var = rawA(:,str2num(assmpt_info(m).row));
        pos = strcat(char(65+str2num(assmpt_info(m).row)),'3');
        xlswrite(output_file,var,assmpt_info(m).co,pos);
    end

end

end


end

