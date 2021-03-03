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
vars.def = 'Foreign demand (GDP growth)';
varsT = [varsT,vars];

%%
%Expenditure side - annual percentage changes
vars.nameq = 'GC';
vars.name = {'GCA'};
vars.def = 'Private consumption expenditure';
varsT = [varsT,vars];

vars.nameq = {'GCG'};
vars.name = {'GCGA'};
vars.def = 'Government consumption expenditure';
varsT = [varsT,vars];

% vars.name = {'GI+GIG'};
% vars.def = 'Gross fixed capital formation';
% varsT = [varsT,vars];

vars.name = {'GIGA'};
vars.def = 'Gross fixed capital formation - General Government';
varsT = [varsT,vars];

vars.name = {'GIA'};
vars.def = 'Gross fixed capital formation - Private sector';
varsT = [varsT,vars];
% 
% vars.name = {''};
% vars.def = 'Domestic demand (incl. Inventories)';
% varsT = [varsT,vars];

vars.name = {'GXA'};
vars.def = 'Exports of goods and services';
varsT = [varsT,vars];

vars.name = {'GMTOTA'};
vars.def = 'Imports of goods and services';
varsT = [varsT,vars];

vars.name = {'GYOBSA'};
vars.def = 'Gross domestic product at market prices';
varsT = [varsT,vars];
%%
% Expenditure side deflators - annual percentage changes
vars.name = {'PHICVATA'};
vars.def = 'Private consumption expenditure';
varsT = [varsT,vars];

vars.name = {'PHIGA'};
vars.def = 'Government consumption expenditure';
varsT = [varsT,vars];

vars.name = {'PHIIGA'};
vars.def = 'Gross fixed capital formation - General Government';
varsT = [varsT,vars];

vars.name = {'PHIIA'};
vars.def = 'Gross fixed capital formation - Private sector';
varsT = [varsT,vars];

vars.name = {'PHIXA'};
vars.def = 'Exports of good and services';
varsT = [varsT,vars];

vars.name = {'PHIMTOTA'};
vars.def = 'Imports of good and services';
varsT = [varsT,vars];

vars.name = {'PHIYOBSA'};
vars.def = 'Gross domestic product at market prices';
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
vars.def = 'External balance of goods and services';
varsT = [varsT,vars];

%%
%General Government account - % of GDP

vars.name = {'RGYA'};
vars.def = 'Total government revenue';
varsT = [varsT,vars];

% vars.name = {'GYA'};
% vars.def = 'Total government expenditure';
% varsT = [varsT,vars];
% 
% vars.name = {''}
% vars.def = 'Net lending (+)/ net borrowing (-)'
% varsT = [varsT,vars]

vars.name = {'BGYA'};
vars.def = 'General government consolidated gross debt';
varsT = [varsT,vars];

vars.name = {'DFETOTA'};
vars.def = 'Total government discretionary fiscal effort';
varsT = [varsT,vars];

vars.name = {'GTA'};
vars.def = 'Government Transfers';
varsT = [varsT,vars];
%%
%Disposable Income; Saving Rate
vars.name = {'GDIA'};
vars.def = 'Nominal Disposable Income';
varsT = [varsT,vars];

vars.name = {'SRA'};
vars.def = 'Saving rate ';
varsT = [varsT,vars];
