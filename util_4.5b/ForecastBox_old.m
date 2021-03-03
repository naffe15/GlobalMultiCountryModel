
function ForecastBox(path, Datafile, year, M_, ncfrcst, mydate, annual)
% This function prepares the shock decomposition box for DG EcFin
% Input:
%   path: path of the experiment
%   Datafile: data file containing the forecasts (including its path)
%   year: forecasted year 
%
% Output:
%   ForecastBox.xls file in the current folder
%
% Example: from the GM3 folder
% ForecastBox('.\2017_03_25_Aytrends_levels_inomss_est_bw4_abw_est_muy_phi0_y0_snap_Bquat_est', '.\2017_03_25_Aytrends_levels_inomss_est_bw4_abw_est_muy_phi0_y0_snap_Bquat_est\src\GM2016Q3.xlsx', 2017)
    
%    year = 2017;
%    Datafile = 'GM2016Q3.xlsx';
%    path = 'H:\IAZZCLUSTER\Global_Estimated_Model\GM\GM3\2017_03_25_Aytrends_levels_inomss_est_bw4_abw_est_muy_phi0_y0_snap_Bquat_est';

if nargin<5 || isempty(ncfrcst),
    ncfrcst=8;
end
if nargin<6 || isempty(mydate),
    mydate=dates('2016Q4');
end
if nargin<7,
    annual=0;
end
idir = findstr(filesep,path);
configu = path(idir(end-1)+1:idir(end)-1);
country = configu(end-1:end);
if ~any(strcmp(country,{'DE','ES','IT','FR'}))
    country = 'EA';
    [Data, Headers, Raw] = xlsread(sprintf('%s', Datafile), 'EA19');
else
    [Data, Headers, Raw] = xlsread(sprintf('%s', Datafile), country);
end
%    country = 'NULL';
%    for i=1:length(countries)    
%        if( strcmp(countries{i}, 'DE') || strcmp(countries{i}, 'FR') || strcmp(countries{i}, 'IT') || strcmp(countries{i}, 'ES') || strcmp(countries{i}, 'NL') )
%            country = countries{i};
%        end
%    end
%    if strcmp(country, 'NULL')
%        country = 'EA';
%        [Data Headers] = xlsread(sprintf('%s', Datafile), 'EA19');
%    else
%        [Data Headers] = xlsread(sprintf('%s', Datafile), country);
%    end

    % Get long-run trend  
    if nargin==3
    disp('Loading gemc_results.mat ...')
    load(sprintf('%s\\gemc_results.mat', path),'M_');
    end
    if annual
        LongTermTrend = 100*get_param_by_name('GYTREND0');
    else
        LongTermTrend = 100*get_param_by_name('GYTREND0')*4;
    end
    
    Cols = 1:size(Data,2);
    try
        indx = strmatch(int2str(year),Headers(1,:)');
        colNo = Cols(indx(1));
        GDP_Real_EA = Raw{169,colNo};
        if ischar(GDP_Real_EA),
            GDP_Real_EA=NaN;
        end
        GDP_Real_EA_1 = Raw{169,colNo-1};
    catch
        indx = find(Data(1,:)==year);
        colNo = Cols(indx(1));
        GDP_Real_EA = Data(169,colNo);
        GDP_Real_EA_1 = Data(169,colNo-1);
    end
    % GDP_Real_EA(isnan(GDP_Real_EA))=[];
    EcFinForecast = 100*(GDP_Real_EA-GDP_Real_EA_1)/GDP_Real_EA_1;

    % Select input files
    if annual
        Forecast_File = sprintf('%s\\%s', path, ['gemc_shock_decomposition 2-step ahead forecast (given 2016Y) group ' country '.xls']);
        CondForecast_File = sprintf('%s\\%s', path, ['gemc_shock_decomposition 2-step ahead conditional forecast (given 2016Y) group ' country '.xls']);
        [F_Data, F_Headers] = xlsread(Forecast_File,['GYOBS_' country]);
        [CF_Data, CF_Headers] = xlsread(CondForecast_File,['GYOBS_' country] );
    else
        Forecast_File = sprintf('%s\\%s', path, ['gemc_shock_decomposition ' int2str(ncfrcst) '-step ahead forecast (given ' char(mydate) ') group ' country '.xls']);
        CondForecast_File = sprintf('%s\\%s', path, ['gemc_shock_decomposition ' int2str(ncfrcst) '-step ahead conditional forecast (given ' char(mydate) ') group ' country '.xls']);
        [F_Data, F_Headers] = xlsread(Forecast_File,['GYOBSA_' country]);
        [CF_Data, CF_Headers] = xlsread(CondForecast_File,['GYOBSA_' country] );
    end
    
    exassmptyes = dir('*assmpt*.xls');
    if ~isempty(exassmptyes)
        Forecast_File_assmpt = sprintf('%s\\%s', path, ['gemc_shock_decomposition assmpt realtime (rolling) group ' country '.xls']);
        if annual
            CondForecast_File_assmpt = sprintf('%s\\%s', path, ['gemc_shock_decomposition assmpt 2-step ahead conditional forecast (given 2016Y) group ' country '.xls']);
            [Fassmpt_Data, Fassmpt_Headers] = xlsread(Forecast_File_assmpt,['GYOBS_' country]);
            [CFassmpt_Data, CFassmpt_Headers] = xlsread(CondForecast_File_assmpt,['GYOBS_' country] );
        else
            CondForecast_File_assmpt = sprintf('%s\\%s', path, ['gemc_shock_decomposition assmpt ' int2str(ncfrcst) '-step ahead conditional forecast (given ' char(mydate) ') group ' country '.xls']);
            [Fassmpt_Data, Fassmpt_Headers] = xlsread(Forecast_File_assmpt,['GYOBSA_' country]);
            [CFassmpt_Data, CFassmpt_Headers] = xlsread(CondForecast_File_assmpt,['GYOBSA_' country] );
        end
    end
    
    
    % Load the data
    % Line 1 = 2016, 2 = 2017
    % Column 1=year, 2=TFP EA, 3=Fiscal EA, 4=Monetary EA, 5=Price Mark-up EA, 6=Bond premium EA vs RoW, 7=Bond premia others vs RoW, 8=Private savings shock EA, 
    %   9=Investment risk premium EA, 10=Wage Mark-up EA, 11=Labor demand shock EA, 12=Other shocks EA, 13=Trade shocks, 14=Shocks US, 15=Shocks RoW, 16=Oil, 
    %   17=Flight to Safety, 18=Others, 19=Init Smoot Var
    Rows = 1:size(F_Data,1);
    rowNo = Cols(F_Data(:,1)==year);
    % Table computations
    TFP         = [ F_Data(rowNo,2) CF_Data(rowNo,2) ];
    indx = [strmatch('Price Mark-up',F_Headers(1,:)) strmatch('Wage Mark-up',F_Headers(1,:)) strmatch('Labor cost',F_Headers(1,:))];
    
    LabGoodMkt  = [ sum(F_Data(rowNo,indx)) sum(CF_Data(rowNo,indx)) ]; % Price Mark-up + Wage Mark-up + Labor demand shock
   
    indx = [strmatch('Oil',F_Headers(1,:)) ];
    
%     if isempty(indx)
%         indx = [strmatch('Commodities',F_Headers(1,:))];
%     end
    
    Oil         = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ]; 

    indx = [strmatch('Non Oil',F_Headers(1,:)) ];
     if ~isempty(indx)
        NonOil         = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
     end
    indx = [strmatch('Private savings',F_Headers(1,:)) strmatch('Flight to safety',F_Headers(1,:))];
    Consumption = [ sum(F_Data(rowNo,indx)) sum(CF_Data(rowNo,indx)) ]; %Private savings shock + Flight to Safety
    indx = strmatch('Investment',F_Headers(1,:));
    Investment  = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ]; 
    indx = strmatch('Fiscal',F_Headers(1,:));
    Fiscal      = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ]; 
    indx = strmatch('Monetary',F_Headers(1,:));
    Monetary    = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ]; 
    
    
    if strcmpi(configu,'gm3')    
 %   if strcmp(country, 'NULL')
           US  = [ F_Data(rowNo,14) CF_Data(rowNo,14) ]; 
    else
        US = [nan nan];
    end
    
    indx = strmatch('Shocks RoW',F_Headers(1,:));
    RoW         = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
    indx = strmatch('Trade shocks',F_Headers(1,:));
    Trade       = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
    indx = strmatch('Bond premium',F_Headers(1,:));
    FX          = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
    
    if exist('NonOil') 
        Tab_Data_tmp = [ TFP; LabGoodMkt; Oil; NonOil; [NaN NaN]; [NaN NaN]; Consumption; Investment; Fiscal; Monetary; [NaN NaN]; US; RoW; Trade; FX; [NaN NaN]; [NaN NaN] ];
    else
        Tab_Data_tmp = [ TFP; LabGoodMkt; Oil;  [NaN NaN]; [NaN NaN]; Consumption; Investment; Fiscal; Monetary; [NaN NaN]; US; RoW; Trade; FX; [NaN NaN]; [NaN NaN] ];
    end
    
    Tab_Data = nan(size(Tab_Data_tmp,1)+2,3);
    Tab_Data(3:size(Tab_Data,1),1:3) = 100*[Tab_Data_tmp (sum(Tab_Data_tmp'))' ];
    
    Tab_Data(2,3) =LongTermTrend;
    Tab_Data(size(Tab_Data_tmp,1)+1,3) = EcFinForecast - LongTermTrend - sum_nan(Tab_Data(3:size(Tab_Data_tmp,1),3));
    Tab_Data(size(Tab_Data,1),3) = EcFinForecast;
    
    
    if ~isempty(exassmptyes)
        TFPassmpt         = [ Fassmpt_Data(rowNo,2) ];%;CFassmpt_Data(rowNo,2) ];
        indx = [strmatch('Price Mark-up',F_Headers(1,:)) strmatch('Wage Mark-up',F_Headers(1,:)) strmatch('Labor cost',F_Headers(1,:))];
        
        LabGoodMktassmpt  = [ sum(Fassmpt_Data(rowNo,indx))];% sum(CF_Data(rowNo,indx)) ]; % Price Mark-up + Wage Mark-up + Labor demand shock
        indx = strmatch('Oil',F_Headers(1,:));
%         if isempty(indx)
%             indx = [strmatch('Commodities',F_Headers(1,:))];
%         end
        Oilassmpt         = [ Fassmpt_Data(rowNo,indx)];% CFassmpt_Data(rowNo,indx) ];
        indx = [strmatch('Non Oil',F_Headers(1,:)) ];
        if ~isempty(indx)
            NonOilassmpt         = [ Fassmpt_Data(rowNo,indx)  ];
        end
        indx = [strmatch('Private savings',F_Headers(1,:)) strmatch('Flight to safety',F_Headers(1,:))];
        Consumptionassmpt = [ sum(Fassmpt_Data(rowNo,indx))];% sum(CFassmpt_Data(rowNo,indx)) ]; %Private savings shock + Flight to Safety
        indx = strmatch('Investment',F_Headers(1,:));
        Investmentassmpt  = [ Fassmpt_Data(rowNo,indx)];% CFassmpt_Data(rowNo,indx) ];
        indx = strmatch('Fiscal',F_Headers(1,:));
        Fiscalassmpt      = [ Fassmpt_Data(rowNo,indx)];% CFassmpt_Data(rowNo,indx) ];
        indx = strmatch('Monetary',F_Headers(1,:));
        Monetaryassmpt    = [ Fassmpt_Data(rowNo,indx)];% CFassmpt_Data(rowNo,indx) ];
        indx = strmatch('Smoot Var',F_Headers(1,:));
        SmoothVarAssmpt = [ Fassmpt_Data(rowNo,indx) ];
        RealGDPAssmpt = SmoothVarAssmpt*100 + LongTermTrend;
        if strcmpi(configu,'gm3')
            %   if strcmp(country, 'NULL')
            USassmpt  = [Fassmpt_Data(rowNo,14)];% CFassmpt_Data(rowNo,14) ];
        else
            USassmpt = [nan];% nan];
        end
        
        indx = strmatch('Shocks RoW',F_Headers(1,:));
        RoWassmpt         = [ Fassmpt_Data(rowNo,indx)];% CFassmpt_Data(rowNo,indx) ];
        indx = strmatch('Trade shocks',F_Headers(1,:));
        Tradeassmpt       = [ Fassmpt_Data(rowNo,indx)];% CFassmpt_Data(rowNo,indx) ];
        indx = strmatch('Bond premium',F_Headers(1,:));
        FXassmpt          = [ Fassmpt_Data(rowNo,indx)];% CFassmpt_Data(rowNo,indx) ];
        if exist('NonOil') 
            Tab_Data_tmp_assmpt = [ TFPassmpt; LabGoodMktassmpt; Oilassmpt; NonOilassmpt; [NaN]; [NaN ]; Consumptionassmpt; Investmentassmpt; Fiscalassmpt; Monetaryassmpt; [NaN ]; USassmpt; RoWassmpt; Tradeassmpt; FXassmpt; [NaN ]; [NaN] ];
        else
            Tab_Data_tmp_assmpt = [ TFPassmpt; LabGoodMktassmpt; Oilassmpt;  [NaN]; [NaN ]; Consumptionassmpt; Investmentassmpt; Fiscalassmpt; Monetaryassmpt; [NaN ]; USassmpt; RoWassmpt; Tradeassmpt; FXassmpt; [NaN ]; [NaN] ];
        end
        
        Tab_Data_assmpt = nan(size(Tab_Data_tmp_assmpt,1)+2,3);
        Tab_Data_assmpt(3:size(Tab_Data,1),1) = 100*[Tab_Data_tmp_assmpt];
    
        Tab_Data_assmpt(2,1) =LongTermTrend;
        Tab_Data_assmpt(size(Tab_Data_tmp_assmpt,1)+1,1) = RealGDPAssmpt - LongTermTrend - sum_nan(Tab_Data_assmpt(3:size(Tab_Data_tmp_assmpt,1),1));
        Tab_Data_assmpt(size(Tab_Data_assmpt,1),1) = RealGDPAssmpt;
        
        
           
        Tab_Data_TOT = [Tab_Data(:,1) Tab_Data_assmpt(:,1)-Tab_Data(:,1) Tab_Data_assmpt(:,1) Tab_Data(:,3)-Tab_Data_assmpt(:,1) Tab_Data(:,3)];
    end
    
    
    % Table preparation
    box = cell(15, 4);
    box(1,2) = {['JRC forecast box ' int2str(year)]};
    box(2,2:4) = {'Historical', 'Forecast', 'Total'};
    if exist('NonOil')
        box(3:size(Tab_Data,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Oil', 'Non Oil', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Foreign', 'US', 'RoW', 'Trade', 'Exchange rate', 'Others', 'Real GDP growth (from forecast)'};
    else
        box(3:size(Tab_Data,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Oil', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Foreign', 'US', 'RoW', 'Trade', 'Exchange rate', 'Others', 'Real GDP growth (from forecast)'};
    end
    box(3:size(Tab_Data,1)+2, 2:4) = num2cell(Tab_Data);

     if ~isempty(exassmptyes)
        box_assmpt = cell(15, 6);
        box_assmpt(1,2) = {['JRC forecast box ' int2str(year)]};
        box_assmpt(2,2:6) = {'Historical','Assumptions', 'Total (Historical + Assumptions)', 'Forecast (not assumptions)', 'Total'};
        if exist('NonOil')
            box_assmpt(3:size(Tab_Data_assmpt,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Oil', 'Non Oil', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Foreign', 'US', 'RoW', 'Trade', 'Exchange rate', 'Others', 'Real GDP growth (from forecast)'};
        else
            box_assmpt(3:size(Tab_Data_assmpt,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Oil', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Foreign', 'US', 'RoW', 'Trade', 'Exchange rate', 'Others', 'Real GDP growth (from forecast)'};
        end
        box_assmpt(3:size(Tab_Data_assmpt,1)+2, 2:6) = num2cell(Tab_Data_TOT); 
     end
    
    
    % Export in Excel
    sheet = 'Conditional(No assumption)';
    xlswrite(sprintf('%s\\ForecastBox.xls', path), box, sheet)
if ~isempty(exassmptyes)
    sheet = 'Assmption';
     xlswrite(sprintf('%s\\ForecastBox.xls', path), box_assmpt, sheet);
end