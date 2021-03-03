clear all
cd Z:\Global_Estimated_Model\Tools\Arima_GM
mypath = pwd;

addpath z:\SVN\captain\captain\
addpath Z:\svn\util\
addpath Z:\svn\dyn_util\util_4.5a\
addpath Z:\dynare_git\dynare_4.5_JRC\matlab\
addpath Z:\Global_Estimated_Model\Tools\Arima_GM


%% GM2
cd Z:\Global_Estimated_Model\GM\GM2\2017_11_16_SIGMAZ_AR2_GPC0_rel_price_3ndstorage_noARMAINOM_TRestim
load gemc_results 
%load pdata GPOP0 GYTREND0
gpop_ea = oo_.SmoothedVariables.GPOP_EA - get_mean('GPOP_EA');
gactr_ea = oo_.SmoothedVariables.GACTR_EA - get_mean('GACTR_EA');
gaytrend_ea = oo_.SmoothedVariables.GAYTREND_EA - get_mean('GAYTREND_EA');
gpop_rw = oo_.SmoothedVariables.GPOP_RoW - get_mean('GPOP_RoW');
gaytrend_rw = oo_.SmoothedVariables.GAYTREND_RoW - get_mean('GAYTREND_RoW');

%% GERMANY
cd Z:\Global_Estimated_Model\GM\GM3-DE\2017_11_02_AF17_3rd_GP0_FWDtradeDL
load gemc_results 
%load pdata GPOP0 GYTREND0
gpop_de = oo_.SmoothedVariables.GPOP_DE - get_mean('GPOP_DE');
gactr_de = oo_.SmoothedVariables.GACTR_DE - get_mean('GACTR_DE');
gaytrend_de = oo_.SmoothedVariables.GAYTREND_DE - get_mean('GAYTREND_DE');

%% FRANCE
cd Z:\Global_Estimated_Model\GM\GM3-FR\2017_11_02_AF17_3rd_GP0_FWD_SIGMAZ_tradeDL
load gemc_results 
%load pdata GPOP0 GYTREND0
gpop_fr = oo_.SmoothedVariables.GPOP_FR - get_mean('GPOP_FR');
gactr_fr = oo_.SmoothedVariables.GACTR_FR - get_mean('GACTR_FR');
gaytrend_fr = oo_.SmoothedVariables.GAYTREND_FR - get_mean('GAYTREND_FR');

%% ITALY
cd Z:\Global_Estimated_Model\GM\GM3-IT\2017_11_03_AF17_3rd_GP0_FWD_SIGMAZ_tradeDL_test
load gemc_results
%load pdata GPOP0 GYTREND0
gpop_it = oo_.SmoothedVariables.GPOP_IT - get_mean('GPOP_IT');
gactr_it = oo_.SmoothedVariables.GACTR_IT - get_mean('GACTR_IT');
gaytrend_it = oo_.SmoothedVariables.GAYTREND_IT - get_mean('GAYTREND_IT');

%% SPAIN
cd Z:\Global_Estimated_Model\GM\GM3-ES\2017_11_02_AF17_3rd_GP0_FWD_SIGMAZ_tradeDL
load gemc_results 
%load pdata GPOP0 GYTREND0
gpop_es = oo_.SmoothedVariables.GPOP_ES - get_mean('GPOP_ES');
gactr_es = oo_.SmoothedVariables.GACTR_ES - get_mean('GACTR_ES');
gaytrend_es = oo_.SmoothedVariables.GAYTREND_ES - get_mean('GAYTREND_ES');

Vars         = [gpop_ea gactr_ea gaytrend_ea gpop_rw gaytrend_rw gpop_de gactr_de gaytrend_de gpop_fr gactr_fr gaytrend_fr gpop_it gactr_it gaytrend_it gpop_es gactr_es gaytrend_es];
VarsNames    = {'gpop_ea', 'gactr_ea', 'gaytrend_ea', 'gpop_rw', 'gaytrend_rw', 'gpop_de', 'gactr_de', 'gaytrend_de', 'gpop_fr', 'gactr_fr', 'gaytrend_fr', 'gpop_it', 'gactr_it', 'gaytrend_it', 'gpop_es', 'gactr_es', 'gaytrend_es'};

%% Load the results on the ARIMA(p,0,q) models
    for i = 1 : size(Vars,2)

      [p,q,param_ar,param_ma,se,res,irf] = best_arma(Vars(:,i));
      Results.(VarsNames{i}).p = p;
      Results.(VarsNames{i}).q = q;
      Results.(VarsNames{i}).param_ar = param_ar;
      Results.(VarsNames{i}).param_ma = param_ma;
      Results.(VarsNames{i}).se = se;
      Results.(VarsNames{i}).res = res;
      Results.(VarsNames{i}).irf = irf;

    end

%% Plot the IRFs of the ARIMA(p,0,q) models
    for i = 1 : size(Vars,2)
        
    figure,
    subplot(211),plot(Results.(VarsNames{i}).irf)
    formatSpec = 'IRF of %s %s';
    A1 = VarsNames{i}(1:end-3);
    A2 = VarsNames{i}(end-1:end);
    str = sprintf(formatSpec,A1,A2)
    title(str)
    subplot(212),plot(cumsum(Results.(VarsNames{i}).irf))
    formatSpec = 'Cumulated IRF of %s %s';
    A1 = VarsNames{i}(1:end-3);
    A2 = VarsNames{i}(end-1:end);
    str = sprintf(formatSpec,A1,A2)
    title(str)
       
    end
    
   
%% Save the results in an Excel file
cd(mypath)
names={'AR order','MA order','AR parameters','MA parameters','SE of residuals','Residuals','IRF'};

for i = 1 : size(Vars,2)
    if ~ismac
        [STATUS,MESSAGE] = xlswrite('Arima_GM.xlsx',names,VarsNames{i},'A1:G1')
        [STATUS,MESSAGE] = xlswrite('Arima_GM.xlsx',Results.(VarsNames{i}).p,VarsNames{i},'A2')
        [STATUS,MESSAGE] = xlswrite('Arima_GM.xlsx',Results.(VarsNames{i}).q,VarsNames{i},'B2')
        [STATUS,MESSAGE] = xlswrite('Arima_GM.xlsx',Results.(VarsNames{i}).param_ar',VarsNames{i},'C2')
        [STATUS,MESSAGE] = xlswrite('Arima_GM.xlsx',Results.(VarsNames{i}).param_ma',VarsNames{i},'D2')
        [STATUS,MESSAGE] = xlswrite('Arima_GM.xlsx',Results.(VarsNames{i}).se,VarsNames{i},'E2')
        [STATUS,MESSAGE] = xlswrite('Arima_GM.xlsx',Results.(VarsNames{i}).res,VarsNames{i},'F2')
        [STATUS,MESSAGE] = xlswrite('Arima_GM.xlsx',Results.(VarsNames{i}).irf,VarsNames{i},'G2')
    else
        [STATUS,MESSAGE] = xlswrite_MACOS('Arima_GM.xlsx',names,VarsNames{i},'A1:G1')
        [STATUS,MESSAGE] = xlswrite_MACOS('Arima_GM.xlsx',Results.(VarsNames{i}).p,VarsNames{i},'A2')
        [STATUS,MESSAGE] = xlswrite_MACOS('Arima_GM.xlsx',Results.(VarsNames{i}).q,VarsNames{i},'B2')
        [STATUS,MESSAGE] = xlswrite_MACOS('Arima_GM.xlsx',Results.(VarsNames{i}).param_ar',VarsNames{i},'C2')
        [STATUS,MESSAGE] = xlswrite_MACOS('Arima_GM.xlsx',Results.(VarsNames{i}).param_ma',VarsNames{i},'D2')
        [STATUS,MESSAGE] = xlswrite_MACOS('Arima_GM.xlsx',Results.(VarsNames{i}).se,VarsNames{i},'E2')
        [STATUS,MESSAGE] = xlswrite_MACOS('Arima_GM.xlsx',Results.(VarsNames{i}).res,VarsNames{i},'F2')
        [STATUS,MESSAGE] = xlswrite_MACOS('Arima_GM.xlsx',Results.(VarsNames{i}).irf,VarsNames{i},'G2')
    end

end

