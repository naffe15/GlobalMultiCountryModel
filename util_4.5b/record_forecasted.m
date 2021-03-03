function record_forecasted(nfrcst,M_,oo_, CountryName)

if nargin <4
    CountryName='';
end


outFV=ForecastedVariablesFn(nfrcst, M_,oo_,1);


Forecasted=outFV.ForecastedVariables;

% function record_housing_script(xlsname, Sigma0,Params0)
% usage record_housing_script('record_housing_EU.xls', Sigma0,Params0)
% global M_

% xlsname = strcat('record_forecasted_', CountryName, '.xls');
% % M_.Sigma_e=Sigma0;
% % M_.params=Params0;
% % print_params_shocks_xls(xlsname,'All parameters & shocks');
% % print_rmses_xls(rmse,r2,lnam,rmseK,lnamK,['..\',xlsname],'RMSE''s');

print_forecast_xls(Forecasted, M_, oo_,CountryName)

% print_forecast_xls(Forecasted,fname,sheet,'Forecasted');