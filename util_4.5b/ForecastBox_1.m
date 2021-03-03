
function ForecastBox_vars(path, Datafile, year, M_, ncfrcst, mydate, annual,varbox,q2avec,country)
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
if nargin<7 || isempty(annual)
    annual=0;
end
if nargin<10,

idir = findstr(filesep,path);
configu = path(idir(end-1)+1:idir(end)-1);
country = configu(end-1:end);
gm = configu(1:3);
if ~any(strcmp(country,{'DE','ES','IT','FR','NL'}))
    country = 'EA';
end
end

if nargin<8 || isempty(varbox) || isempty(q2avec)
    varbox = 'LYOBS';
    q2avec.qname=['LYOBS_' country];
    q2avec.plot = 1;
    q2avec.gname=['GYOBSA_' country ];
    q2avec.frcst_name='Real GDP';
    
    if nargin==3
        disp('Loading gemc_results.mat ...')
        load(sprintf('%s\\gemc_results.mat', path),'M_');
    end
    q2avec.GYTREND0 = get_param_by_name('GYTREND0');
end

indx =  strmatch([varbox '_' country], {q2avec.qname});
if q2avec(indx).plot==1
    var_box=q2avec(indx).gname;
    LongTermTrend = 100*q2avec(indx).GYTREND0*4;
    if var_box(1:3)=='PHI'
        Frcstvar= q2avec(indx).frcst_name;
    else
        Frcstvar= [q2avec(indx).frcst_name ' growth'];
    end
else
    LongTermTrend = 0;
    var_box=q2avec(indx).name;
    Frcstvar= q2avec(indx).frcst_name;
end


Forecast_File = sprintf('%s\\%s', path, ['gemc_shock_decomposition ' int2str(ncfrcst) '-step ahead forecast (given ' char(mydate) ') group ' country '.xls']);
CondForecast_File = sprintf('%s\\%s', path, ['gemc_shock_decomposition ' int2str(ncfrcst) '-step ahead conditional forecast (given ' char(mydate) ') group ' country '.xls']);
[F_Data, F_Headers] = xlsread(Forecast_File,var_box);
[CF_Data, CF_Headers] = xlsread(CondForecast_File,var_box );


exassmptyes = dir('*assmpt*.xls');
if ~isempty(exassmptyes)
    Forecast_File_assmpt = sprintf('%s\\%s', path, ['gemc_shock_decomposition assmpt realtime (rolling) group ' country '.xls']);
    CondForecast_File_assmpt = sprintf('%s\\%s', path, ['gemc_shock_decomposition assmpt ' int2str(ncfrcst) '-step ahead conditional forecast (given ' char(mydate) ') group ' country '.xls']);
    [Fassmpt_Data, Fassmpt_Headers] = xlsread(Forecast_File_assmpt,var_box);
    [CFassmpt_Data, CFassmpt_Headers] = xlsread(CondForecast_File_assmpt,var_box );
end

Ecfin_Forecast_File = sprintf('%s\\%s', path, ['gemc_annual_frcst_ecfin.xls']);
[Ecfin_F_Data, Ecfin_F_Headers] = xlsread(Ecfin_Forecast_File,var_box);
%     EcFinForecast = 100*(GDP_Real_EA-GDP_Real_EA_1)/GDP_Real_EA_1;
EcFinForecast = 100*Ecfin_F_Data(strmatch('EcFin',Ecfin_F_Headers),find(Ecfin_F_Data(1,:)==year));

% Load the data
% Line 1 = 2016, 2 = 2017
% Column 1=year, 2=TFP EA, 3=Fiscal EA, 4=Monetary EA, 5=Price Mark-up EA, 6=Bond premium EA vs RoW, 7=Bond premia others vs RoW, 8=Private savings shock EA,
%   9=Investment risk premium EA, 10=Wage Mark-up EA, 11=Labor demand shock EA, 12=Other shocks EA, 13=Trade shocks, 14=Shocks US, 15=Shocks RoW, 16=Oil,
%   17=Flight to Safety, 18=Others, 19=Init Smoot Var

%DETAILED
% 1=TFP 2=Fiscal 3=Monetary 4=Price Mark-up 5=Bond premium RoW vs EA
% 6=Private savings shock EA 7=Investment risk premium EA   8=Wage Mark-up E
% 9= Labor cost shock EA    10=Other shocks EA  11=Shocks RoW   12=Trade sho
% 13=Oil    14=Non-oil 14= Flight to safety     15= Others + Initial Values
% 16=smoothed

Rows = 1:size(F_Data,1);
rowNo = find(F_Data(:,1)==year);
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
indx = [strmatch('Commodities demand',F_Headers(1,:)) ];
if ~isempty(indx)
    Commodities_demand         = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
else
    Commodities_demand         = [nan nan];
end


indx = [strmatch('Private savings',F_Headers(1,:)) strmatch('Flight to safety',F_Headers(1,:))];
Consumption = [ sum(F_Data(rowNo,indx)) sum(CF_Data(rowNo,indx)) ]; %Private savings shock + Flight to Safety
indx = strmatch('Investment',F_Headers(1,:));
Investment  = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
indx = strmatch('Fiscal',F_Headers(1,:));
Fiscal      = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
indx = strmatch('Monetary',F_Headers(1,:));
Monetary    = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];

gm = 'gm-emu';
if ~isempty(evalin('base','strmatch(''Y_US'',M_.endo_names)'))
    %   if strcmp(country, 'NULL')
    US  = [ F_Data(rowNo,14) CF_Data(rowNo,14) ];
    gm = 'gm3';
else
    US = [nan nan];
    if isempty(evalin('base','strmatch(''Y_RoW'',M_.endo_names)'))
        gm = 'gm1';
    end
end

if strcmpi(gm,'gm1')==0
indx = strmatch('Shocks RoW',F_Headers(1,:));
RoW         = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
indx = strmatch('Trade shocks',F_Headers(1,:));
Trade       = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
indx = strmatch('Bond premium',F_Headers(1,:));
FX          = [ F_Data(rowNo,indx) CF_Data(rowNo,indx) ];
end

if strcmpi(gm,'gm1')==0
if exist('NonOil')
    Tab_Data_tmp = [ TFP; LabGoodMkt; Oil; NonOil; [NaN NaN]; [NaN NaN]; Consumption; Investment; Fiscal; Monetary; [NaN NaN]; US; RoW; Trade; FX;Commodities_demand; [NaN NaN]; [NaN NaN] ];
else
    Tab_Data_tmp = [ TFP; LabGoodMkt; Oil;  [NaN NaN]; [NaN NaN]; Consumption; Investment; Fiscal; Monetary; [NaN NaN]; US; RoW; Trade; FX;Commodities_demand; [NaN NaN]; [NaN NaN] ];
end
end
Tab_Data = nan(size(Tab_Data_tmp,1)+2,3);
Tab_Data(3:size(Tab_Data,1),1:3) = 100*[Tab_Data_tmp (sum(Tab_Data_tmp'))' ];

Tab_Data_tmp2 = [ F_Data(rowNo,2:end); CF_Data(rowNo,2:end) ]';
Tab_Data2= nan(size(Tab_Data_tmp2,1)+2,3);
Tab_Data2(3:size(Tab_Data2,1),1:3) = 100*[Tab_Data_tmp2 (sum(Tab_Data_tmp2'))' ];

%others
Tab_Data(2,3) =LongTermTrend;
Tab_Data(size(Tab_Data_tmp,1)+1,3) = EcFinForecast - LongTermTrend - sum_nan(Tab_Data(3:size(Tab_Data_tmp,1),3));
Tab_Data(size(Tab_Data,1),3) = EcFinForecast;

Tab_Data2(size(Tab_Data_tmp2,1)+1,3) = EcFinForecast - LongTermTrend - sum_nan(Tab_Data2(3:size(Tab_Data_tmp2,1),3));


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
    indx = [strmatch('Commodities demand',F_Headers(1,:)) ];
    if ~isempty(indx)
        Commodities_demandassmpt         = [ Fassmpt_Data(rowNo,indx)  ];
    else
        Commodities_demandassmpt         = [ nan ];
       
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
    if strcmpi(gm,'gm3')
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
        Tab_Data_tmp_assmpt = [ TFPassmpt; LabGoodMktassmpt; Oilassmpt; NonOilassmpt; [NaN]; [NaN ]; Consumptionassmpt; Investmentassmpt; Fiscalassmpt; Monetaryassmpt; [NaN ]; USassmpt; RoWassmpt; Tradeassmpt; FXassmpt;Commodities_demandassmpt; [NaN ]; [NaN] ];
    else
        Tab_Data_tmp_assmpt = [ TFPassmpt; LabGoodMktassmpt; Oilassmpt;  [NaN]; [NaN ]; Consumptionassmpt; Investmentassmpt; Fiscalassmpt; Monetaryassmpt; [NaN ]; USassmpt; RoWassmpt; Tradeassmpt; FXassmpt;Commodities_demandassmpt; [NaN ]; [NaN] ];
    end
    
    Tab_Data_assmpt = nan(size(Tab_Data_tmp_assmpt,1)+2,3);
    Tab_Data_assmpt(3:size(Tab_Data,1),1) = 100*[Tab_Data_tmp_assmpt];
    
    Tab_Data_assmpt(2,1) =LongTermTrend;
    Tab_Data_assmpt(size(Tab_Data_tmp_assmpt,1)+1,1) = RealGDPAssmpt - LongTermTrend - sum_nan(Tab_Data_assmpt(3:size(Tab_Data_tmp_assmpt,1),1));
    Tab_Data_assmpt(size(Tab_Data_assmpt,1),1) = RealGDPAssmpt;
    
    Tab_Data_tmp_assmpt2 = [ Fassmpt_Data(rowNo,2:end) ];
    Tab_Data_assmpt2 = nan(size(Tab_Data_tmp_assmpt2,1)+2,3);
    Tab_Data_assmpt2(3:size(Tab_Data2,1),1) = 100*[Tab_Data_tmp_assmpt2];
    
    Tab_Data_assmpt2(size(Tab_Data_tmp_assmpt2,1)+1,1) = RealGDPAssmpt - LongTermTrend - sum_nan(Tab_Data_assmpt2(3:size(Tab_Data_tmp_assmpt2,1),1));
    Tab_Data_assmpt2(size(Tab_Data_assmpt2,1),1) = RealGDPAssmpt;
    
    
    Tab_Data_TOT = [Tab_Data(:,1) Tab_Data_assmpt(:,1)-Tab_Data(:,1) Tab_Data_assmpt(:,1) Tab_Data(:,3)-Tab_Data_assmpt(:,1) Tab_Data(:,3)];
    Tab_Data_TOT(end-1:end,4)=nan;
    
    Tab_Data_TOT2 = [Tab_Data2(:,1) Tab_Data_assmpt2(:,1)-Tab_Data2(:,1) Tab_Data_assmpt2(:,1) Tab_Data2(:,3)-Tab_Data_assmpt2(:,1) Tab_Data2(:,3)];
    Tab_Data_TOT(end-1:end,4)=nan;
    
end


% Table preparation
box = cell(15, 4);
box(1,2) = {['JRC forecast box ' int2str(year)]};
box(2,2:4) = {'Historical', 'Forecast', 'Total'};

if strcmpi(gm,'gm1')==0
if exist('NonOil')
    box(3:size(Tab_Data,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Oil', 'Non Oil', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Foreign', 'US', 'RoW', 'Trade', 'Exchange rate', 'Commodities demand','Others', Frcstvar};
else
    box(3:size(Tab_Data,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Oil', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Foreign', 'US', 'RoW', 'Trade', 'Exchange rate', 'Commodities demand','Others', Frcstvar};
end
else
    box(3:size(Tab_Data,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Others', Frcstvar};
end
box(3:size(Tab_Data,1)+2, 2:4) = num2cell(Tab_Data);


box2=cell(19,4);
box2(1,2) = {'Detailed '};
box2(2,2:4) = {'Historical', 'Forecast', 'Total'};
box2(3:size(Tab_Data2,1), 1) = F_Headers(1,2:end);
box2(3:size(Tab_Data2,1), 2:4) = num2cell(Tab_Data2(3:end,:));

if ~isempty(exassmptyes)
    box_assmpt = cell(15, 6);
    box_assmpt(1,2) = {['JRC forecast box ' int2str(year)]};
    box_assmpt(2,2:6) = {'Historical','Assumptions', 'Total (Historical + Assumptions)', 'Forecast (not assumptions)', 'Total'};
    if exist('NonOil')
        box_assmpt(3:size(Tab_Data_assmpt,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Oil', 'Non Oil', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Foreign', 'US', 'RoW', 'Trade', 'Exchange rate', 'Commodities demand','Others', Frcstvar};
    else
        box_assmpt(3:size(Tab_Data_assmpt,1)+2, 1) = {'Supply', 'Long-run trend', 'TFP', 'Labor and good market', 'Oil', 'Demand', 'Domestic', 'Consumption', 'Investment', 'Fiscal', 'Monetary policy', 'Foreign', 'US', 'RoW', 'Trade', 'Exchange rate','Commodities demand', 'Others', Frcstvar};
    end
    box_assmpt(3:size(Tab_Data_assmpt,1)+2, 2:6) = num2cell(Tab_Data_TOT);
    
    box_assmpt2 = cell(19, 6);
    box_assmpt2(1,2) = {['Detailed']};
    box_assmpt2(2,2:6) = {'Historical','Assumptions', 'Total (Historical + Assumptions)', 'Forecast (not assumptions)', 'Total'};
    box_assmpt2(3:size(Tab_Data_assmpt2,1), 1) = F_Headers(1,2:end);
    box_assmpt2(3:size(Tab_Data_assmpt2,1), 2:6) = num2cell(Tab_Data_TOT2(3:end,:));
    
end


% Export in Excel
sheet = [var_box '_Cond(No assumption)'];
xlswrite(sprintf('%s\\ForecastBox.xls', path), box, sheet)
xlswrite(sprintf('%s\\ForecastBox.xls', path), box2(1:end,1:4), sheet,'A24')
xlswrite(sprintf('%s\\ForecastBox.xls', path),    F_Headers(6:end,:), sheet,'A50')


if ~isempty(exassmptyes)
    sheet = [var_box '_Assmption'];
    xlswrite(sprintf('%s\\ForecastBox.xls', path), box_assmpt, sheet);
    xlswrite(sprintf('%s\\ForecastBox.xls', path), box_assmpt2(1:end,1:6), sheet,'A24')
     xlswrite(sprintf('%s\\ForecastBox.xls', path),    F_Headers(6:end,:), sheet,'A50')

    
end