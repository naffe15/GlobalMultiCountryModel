% function record_housing_script(xlsname, Sigma0,Params0)
% usage record_housing_script('record_housing_EU.xls', Sigma0,Params0)
% global M_

xlsname = 'record_GM3.xls';
% M_.Sigma_e=Sigma0;
% M_.params=Params0;
print_params_shocks_xls(['..\',xlsname],'All parameters & shocks');
print_rmses_xls(rmse,r2,lnam,rmseK,lnamK,['..\',xlsname],'RMSE''s');
