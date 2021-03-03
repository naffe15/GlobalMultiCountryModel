% script to make scenarios w.r.t. baseline using ANNUAL assumptions
num = xlsread('src/ea_trade_assmpt.xlsx');
o.T = num(1,:);
o.GMTOT_EA=num(2,:)/100;
o.GX_EA=num(3,:)/100;
f(1)=1;
for k=2:4
    f(k)=f(k-1)*(o.GMTOT_EA(k)+1);
end
b.MTOT_EA = fill_data(NaN(12,1),f(2:end),2,2,[],1,1);
f(1)=1;
for k=2:4
    f(k)=f(k-1)*(o.GX_EA(k)+1);
end
b.X_EA = fill_data(NaN(12,1),f(2:end),2,2,[],1,1);
b.MTOT_EA = fill_data(NaN(16,1),exp(cumsum(o.GMTOT_EA(1:end))),2,2,[],1,1);
b.X_EA = fill_data(NaN(16,1),exp(cumsum(o.GX_EA(1:end))),2,2,[],1,1);
x=[ones; b.X_EA];
b.GX_EA =x(2:end)./x(1:end-1)-1;

constrained_vars_ = [];
constrained_paths_ = zeros(1, 12);
constrained_vars_ = strmatch('GX_EA',M_.endo_names,'exact'); %345;
constrained_paths_ = b.GX_EA'+get_mean('GX_EA') ;
% constrained_paths_(1,1)=tmp2(1);
% constrained_paths_(1,2)=tmp2(2);
% constrained_paths_(1,3)=tmp2(3);
% constrained_paths_(1,4)=tmp2(4);
% constrained_paths_(1,5)=tmp2(5);
% constrained_paths_(1,6)=tmp2(6);
% constrained_paths_(1,7)=tmp2(7);
% constrained_paths_(1,8)=tmp2(8);
shock_name = 'EPS_M_RoW';
% shock_name = 'EPS_UC_RoW';
ind_shock=strmatch(shock_name,M_.exo_names,'exact');
Sigma_e_0=M_.Sigma_e;
M_.Sigma_e = zeros(size(Sigma_e_0,1));
M_.Sigma_e(ind_shock,ind_shock)=Sigma_e_0(ind_shock,ind_shock);
options_cond_fcst_ = struct();
options_cond_fcst_.periods = 12;
options_cond_fcst_.replic = 10;
options_cond_fcst_.parameter_set = 'calibration';
options_cond_fcst_.controlled_varexo = char(shock_name);
imcforecast(constrained_paths_, constrained_vars_, options_cond_fcst_);
load conditional_forecasts forecasts
oo_.jrc.forecast_first_step.Mean = forecasts.cond.Mean;
oo_.jrc.forecast_first_step.ci   = forecasts.cond.ci; 
M_.Sigma_e=Sigma_e_0;
vname = 'LX_EA';
[ya, yass, gya, gyass] = quarterly2annual([forecasts.cond.Mean.(vname)(1)*ones(3,1); forecasts.cond.Mean.(vname)],0,0,2,1,0);
vname = 'LYOBS_EA';
[ya, yass, gya, gyass] = quarterly2annual([forecasts.cond.Mean.(vname)(1)*ones(3,1); forecasts.cond.Mean.(vname)],0,0,2,1,0);
