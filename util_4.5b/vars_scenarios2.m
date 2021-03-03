function varsT=vars_scenarios
varsT=[];
%%
% y,yss,GYTREND0,type,islog,aux
%Technical assumptions

vars.assmpt = 1;
vars.yss = 0;
vars.GYTREND0 = 0;
vars.type = 1;
vars.islog = 0;
vars.aux = 0;
vars.nameq = 'GE_EA';
vars.name = {'GEA_EA'};
vars.def = 'Nominal effective exchange rate(%change)';
varsT = [varsT,vars];

vars.assmpt = 1;
vars.yss = 0;
vars.GYTREND0 = 0;
vars.type = 1;
vars.islog = 0;
vars.aux = 0;
vars.nameq = 'INOM_EA';
vars.name = {'INOMA_EA'};
vars.def = 'Short-term interest rate(%)';
varsT = [varsT,vars];

vars.assmpt = 1;
vars.yss = 1;
vars.GYTREND0 = 0;
vars.type = 2;
vars.islog = 0;
vars.aux = 0;
vars.nameq = 'PHIBRENTOBS_RoW';
vars.name = {'PHIBRENTOBSA_RoW'};
vars.def = 'Oil prices (USD/bbl)';
varsT = [varsT,vars];

vars.assmpt = 1;
vars.yss = get_param_by_name('GYTREND0');
vars.GYTREND0 = 0;
vars.type = 1;
vars.islog = 0;
vars.aux = 0;
vars.nameq = 'GYOBS_RoW';
vars.name = {'GYOBSA_RoW'};
vars.def = 'Growth rate of GDP - RoW';
varsT = [varsT,vars];

%%
%Expenditure side - annual percentage changes
vars.nameq = 'GC';
vars.name = {'GCA'};
vars.def = 'Private Consumption';
varsT = [varsT,vars];

vars.nameq = 'GCS';
vars.name = {'GCSA'};
vars.def = 'CS Private Consumption';
varsT = [varsT,vars];

vars.nameq = 'GCC';
vars.name = {'GCCA'};
vars.def = 'CC Private Consumption';
varsT = [varsT,vars];


vars.nameq = {'GCG'};
vars.name = {'GCGA'};
vars.def = 'Government Consumption';
varsT = [varsT,vars];

% vars.name = {'GI+GIG'};
% vars.def = 'Gross fixed capital formation';
% varsT = [varsT,vars];

vars.name = {'GIGA'};
vars.def = 'Gross Fixed Capital Formation - General Government';
varsT = [varsT,vars];

vars.name = {'GIA'};
vars.def = 'Gross Fixed Capital Formation - Private ';
varsT = [varsT,vars];
% 
% vars.name = {''};
% vars.def = 'Domestic demand (incl. Inventories)';
% varsT = [varsT,vars];

vars.name = {'GXA'};
vars.def = 'Exports of Goods and Services';
varsT = [varsT,vars];

vars.name = {'GMTOTA'};
vars.def = 'Imports of Goods and Services';
varsT = [varsT,vars];

vars.name = {'GYOBSA'};
vars.def = 'Gross Domestic Product at market prices';
varsT = [varsT,vars];

vars.nameq = 'GPYMRKUP';
vars.name = {'GPYMRKUPA'};
vars.def = 'Growth Rate of Price Markup';
varsT = [varsT,vars];
%%
% Expenditure side deflators - annual percentage changes
vars.name = {'PHICVATA'};
vars.def = 'Private Consumption';
varsT = [varsT,vars];

vars.name = {'PHIGA'};
vars.def = 'Government Consumption';
varsT = [varsT,vars];

vars.name = {'PHIIGA'};
vars.def = 'Gross Fixed Capital Formation - General Government';
varsT = [varsT,vars];

vars.name = {'PHIIA'};
vars.def = 'Gross Fixed Capital Formation - Private sector';
varsT = [varsT,vars];

vars.name = {'PHIXA'};
vars.def = 'Exports of Good and Services';
varsT = [varsT,vars];

vars.name = {'PHIMTOTA'};
vars.def = 'Imports of Good and Services';
varsT = [varsT,vars];

vars.name = {'PHIYOBSA'};
vars.def = 'Gross Domestic Product at market prices';
varsT = [varsT,vars];

% vars.name = {''}
% vars.def = 'Terms of trade of goods and services'
% varsT = [varsT,vars]

%%
% Other price indicators

% vars.name = {''}
% vars.def = 'HICP'
% varsT = [varsT,vars]

vars.name = {'PHIWA'};
vars.def = 'Nominal wage inflation';
varsT = [varsT,vars];

vars.name = {'INOMGA'};
vars.def = 'Nominal Interest Rate on Domestic Bonds';
varsT = [varsT,vars];

vars.name = {'INOMGLTA'};
vars.def = 'Long-term gov bonds nominal interest rate';
varsT = [varsT,vars];

vars.name = {'RA'};
vars.def = 'Real interest rate';
varsT = [varsT,vars];

vars.name = {'RSA'};
vars.def = 'Real interest rate on investment';
varsT = [varsT,vars];

vars.name = {'INOMSA'};
vars.def = 'Nominal Interest Rate on Stocks of assets';
varsT = [varsT,vars];


%%
%Supply of goods and services
vars.name = {'CUA'};
vars.def = 'Capacity utilization';
varsT = [varsT,vars];

vars.name = {'FNtNA'};
vars.def = 'Labour hoarding (% total hours)';
varsT = [varsT,vars];

% vars.name = {'TFPA'};
% vars.def = 'TFP';
% varsT = [varsT,vars];

% vars.name = {''}
% vars.def = 'Others (if relevant)'
% varsT = [varsT,vars]

%%
%Labour market indicators

vars.name = {'GNA'};
vars.def = 'Total employment (hours)';
varsT = [varsT,vars];

vars.name = {'YLA'};
vars.def = 'Labor productivity';
varsT = [varsT,vars];

% %%
% %Households+NPISH
% 
% vars.name = {''}
% vars.def = 'Nominal gross disposable income'
% varsT = [varsT,vars]
% 
% vars.name = {''}
% vars.def = 'Real gross disposable income'
% varsT = [varsT,vars]
% 
% vars.name = {''}
% vars.def = 'Saving rate (%)'
% varsT = [varsT,vars]

%%
%External accounts - % of GDP

vars.name = {'TBYA'};
vars.def = 'External balance of goods and services - % of GDP';
varsT = [varsT,vars];


vars.name = {'NFAYA'};
vars.def = 'Net Foreign Asset - % of GDP';
varsT = [varsT,vars];

%%
%General Government account - % of GDP

vars.name = {'RGYA'};
vars.def = 'Total Government Revenues - % of GDP';
varsT = [varsT,vars];

vars.name = {'RGTAUNYA'};
vars.def = '    - Govern. Revenues from labour taxes - % of GDP';
varsT = [varsT,vars];

vars.name = {'RGTAUCYA'};
vars.def = '    - Govern. Revenues from VAT - % of GDP';
varsT = [varsT,vars];

vars.name = {'RGTAUKYA'};
vars.def = '    - Govern. Revenues from corporate taxes - % of GDP ';
varsT = [varsT,vars];

vars.name = {'RGTAXYA'};
vars.def = '    - Govern. Revenues from LS-TAX - % of GDP';
varsT = [varsT,vars];

vars.name = {'RGTAUOILYA'};
vars.def = '    - Govern. Revenues from oil - % of GDP';
varsT = [varsT,vars];

% vars.name = {'GYA'};
% vars.def = 'Total government expenditure';
% varsT = [varsT,vars];
% 
% vars.name = {''}
% vars.def = 'Net lending (+)/ net borrowing (-)'
% varsT = [varsT,vars]

vars.name = {'BGYA'};
vars.def = 'General Government consolidated gross debt - % of GDP';
varsT = [varsT,vars];

vars.name = {'IBYA'};
vars.def = 'General Government interest payments on gross debt - % of GDP';
varsT = [varsT,vars];

vars.name = {'INOMA'};
vars.def = 'General Government Interest Rate on gross debt ';
varsT = [varsT,vars];

vars.name = {'DFETOTA'};
vars.def = 'Total Government Discretionary Fiscal Effort - % of GDP ';
varsT = [varsT,vars];

vars.name = {'PSGYA'};
vars.def = 'General Government Primary Surplus - % of GDP';
varsT = [varsT,vars];

vars.name = {'CGYA'};
vars.def = 'Government Consumption - % of GDP ';
varsT = [varsT,vars];

vars.name = {'GTA'};
vars.def = 'Government Transfers';
varsT = [varsT,vars];

vars.name = {'TYA'};
vars.def = 'Government Transfers - % of GDP ';
varsT = [varsT,vars];

vars.name = {'ZEPS_APC'};
vars.def = 'ZEPS_APC';
varsT = [varsT,vars];

vars.name = {'ZEPS_APG'};
vars.def = 'ZEPS_APG';
varsT = [varsT,vars];

vars.name = {'ZEPS_API'};
vars.def = 'ZEPS_API';
varsT = [varsT,vars];

vars.name = {'ZEPS_APIG'};
vars.def = 'ZEPS_APIG';
varsT = [varsT,vars];

vars.name = {'ZEPS_AY'};
vars.def = 'ZEPS_AY';
varsT = [varsT,vars];

vars.name = {'ZEPS_CU'};
vars.def = 'ZEPS_CU';
varsT = [varsT,vars];

vars.name = {'ZEPS_FQ'};
vars.def = 'ZEPS_FQ';
varsT = [varsT,vars];

vars.name = {'ZEPS_FN'};
vars.def = 'ZEPS_FN';
varsT = [varsT,vars];

vars.name = {'ZEPS_ND'};
vars.def = 'ZEPS_ND';
varsT = [varsT,vars];

vars.name = {'ZEPS_BW'};
vars.def = 'ZEPS_BW';
varsT = [varsT,vars];

vars.name = {'ZEPS_B'};
vars.def = 'ZEPS_B';
varsT = [varsT,vars];

vars.name = {'ZEPS_INOM'};
vars.def = 'ZEPS_INOM';
varsT = [varsT,vars];


vars.name = {'ZEPS_G'};
vars.def = 'ZEPS_G';
varsT = [varsT,vars];

vars.name = {'ZEPS_ACTR'};
vars.def = 'ZEPS_ACTR';
varsT = [varsT,vars];

vars.name = {'ZEPS_HPERE'};
vars.def = 'ZEPS_HPERE';
varsT = [varsT,vars];

vars.name = {'ZEPS_PARTR'};
vars.def = 'ZEPS_PARTR';
varsT = [varsT,vars];

vars.name = {'ZEPS_POP'};
vars.def = 'ZEPS_POP';
varsT = [varsT,vars];

vars.name = {'ZEPS_IG'};
vars.def = 'ZEPS_IG';
varsT = [varsT,vars];


vars.name = {'ZEPS_PX'};
vars.def = 'ZEPS_PX';
varsT = [varsT,vars];

vars.name = {'ZEPS_S'};
vars.def = 'ZEPS_S';
varsT = [varsT,vars];


vars.name = {'ZEPS_T'};
vars.def = 'ZEPS_T';
varsT = [varsT,vars];


vars.name = {'ZEPS_TC'};
vars.def = 'ZEPS_TC';
varsT = [varsT,vars];

vars.name = {'ZEPS_TAUC'};
vars.def = 'ZEPS_TAUC';
varsT = [varsT,vars];

vars.name = {'ZEPS_TAUN'};
vars.def = 'ZEPS_TAUN';
varsT = [varsT,vars];

vars.name = {'ZEPS_TAUK'};
vars.def = 'ZEPS_TAUK';
varsT = [varsT,vars];

vars.name = {'ZEPS_TAX'};
vars.def = 'ZEPS_TAX';
varsT = [varsT,vars];

vars.name = {'ZEPS_U'};
vars.def = 'ZEPS_U';
varsT = [varsT,vars];

vars.name = {'ZEPS_UC'};
vars.def = 'ZEPS_UC';
varsT = [varsT,vars];

vars.name = {'ZEPS_BEN'};
vars.def = 'ZEPS_BEN';
varsT = [varsT,vars];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REA

vars.name = {'GYOBSA_REA'};
vars.def = 'Growth rate of GDP' ;
varsT = [varsT,vars];

vars.name = {'GCA_REA'};
vars.def = 'Growth rate of consumption';
varsT = [varsT,vars];

vars.name = {'PHIYOBSA_REA'};
vars.def = 'Growth rate of GDP price deflator';
varsT = [varsT,vars];

vars.name = {'INOMA_REA'};
vars.def = 'Nominal Interest Rate';
varsT = [varsT,vars];

vars.name = {'TBYA_REA'};
vars.def = 'External balance of goods and services - % of GDP';
varsT = [varsT,vars];


vars.name = {'ZEPS_POP_REA'};
vars.def = 'ZEPS_POP_REA';
varsT = [varsT,vars];

vars.name = {'ZEPS_UC_REA'};
vars.def = 'ZEPS_UC_REA';
varsT = [varsT,vars];

vars.name = {'ZEPS_M_REA'};
vars.def = 'ZEPS_M_REA';
varsT = [varsT,vars];

vars.name = {'ZEPS_Y_REA'};
vars.def = 'ZEPS_Y_REA';
varsT = [varsT,vars];

vars.name = {'ZEPS_PX_REA'};
vars.def = 'ZEPS_PX_REA';
varsT = [varsT,vars];

vars.name = {'ZEPS_A_REA'};
vars.def = 'ZEPS_A_REA';
varsT = [varsT,vars];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ZEPS RoW


% vars.name = {'GYOBSA_RoW'};
% vars.def = 'Growth rate of GDP' ;
% varsT = [varsT,vars];

vars.name = {'GCA_RoW'};
vars.def = 'Growth rate of consumption';
varsT = [varsT,vars];

vars.name = {'PHIYOBSA_RoW'};
vars.def = 'Growth rate of GDP price deflator';
varsT = [varsT,vars];

vars.name = {'INOMA_RoW'};
vars.def = 'Nominal Interest Rate';
varsT = [varsT,vars];

vars.name = {'TBYA_RoW'};
vars.def = 'External balance of goods and services - % of GDP';
varsT = [varsT,vars];



vars.name = {'ZEPS_POP_RoW'};
vars.def = 'ZEPS_POP_RoW';
varsT = [varsT,vars];

vars.name = {'ZEPS_UC_RoW'};
vars.def = 'ZEPS_UC_RoW';
varsT = [varsT,vars];

vars.name = {'ZEPS_M_RoW'};
vars.def = 'ZEPS_M_RoW';
varsT = [varsT,vars];

vars.name = {'ZEPS_Y_RoW'};
vars.def = 'ZEPS_Y_RoW';
varsT = [varsT,vars];

vars.name = {'ZEPS_INOM_RoW'};
vars.def = 'ZEPS_INOM_RoW';
varsT = [varsT,vars];


vars.name = {'ZEPS_A_RoW'};
vars.def = 'ZEPS_A_RoW';
varsT = [varsT,vars];





