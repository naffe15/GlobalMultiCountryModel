function record_script_GM(xlsname, rmse,r2,lnam,rmseK,lnamK)
% usage record_script_GM('record_GM3.xls', rmse,r2,lnam,rmseK,lnamK)
% global M_
if nargin == 0
    disp('record_script_GM(xlsname, rmse,r2,lnam,rmseK,lnamK)')
    return
end

% xlsname = 'record_GM3.xls';
% M_.Sigma_e=Sigma0;
% M_.params=Params0;


%check if record script is getting full
if ~ismac
    [N,~,RAW]=xlsread(['..\',xlsname],'All parameters & shocks');
else
    fix_xlsread_MACOS(['../',xlsname]);
    [N,~,RAW]=xlsread(['../',xlsname],'All parameters & shocks');
end
  if size(N,2)>250
      new_record_script(xlsname);
  end

if ~ismac
    print_params_shocks_xls(['..\',xlsname],'All parameters & shocks');
    print_rmses_xls(rmse,r2,lnam,rmseK,lnamK,['..\',xlsname],'RMSE''s');
else
    print_params_shocks_xls(['../',xlsname],'All parameters & shocks');
    print_rmses_xls(rmse,r2,lnam,rmseK,lnamK,['../',xlsname],'RMSE''s');
end