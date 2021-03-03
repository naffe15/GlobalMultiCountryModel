oo2_ = oo_;
options2_ = options_;
options_.nobs = 89;
options_.shock_decomp=rmfield(options_.shock_decomp,'i_var'); % this is a fix since assmpt stuff is not squeezed!
%%
oo_.realtime_shock_decomposition = oo2_.jrc.assmpt.realtime_shock_decomposition;
oo_.conditional_shock_decomposition = oo2_.jrc.assmpt.conditional_shock_decomposition;
oo_.realtime_conditional_shock_decomposition = oo2_.jrc.assmpt.realtime_conditional_shock_decomposition;
oo_.realtime_forecast_shock_decomposition =oo2_.jrc.assmpt.realtime_forecast_shock_decomposition;
options_ = set_default_plot_shock_decomposition_options(options_);
options_.plot_shock_decomp.detail_plot = 1;
options_.plot_shock_decomp.interactive = 1;
options_.plot_shock_decomp.nodisplay = 1;
options_.plot_shock_decomp.realtime = 1;
options_.plot_shock_decomp.steadystate = 1;
options_.plot_shock_decomp.vintage = 0;
options_.plot_shock_decomp.write_xls = 1;
options_.plot_shock_decomp.fig_name = 'assmpt';
options_.plot_shock_decomp.type = 'aoa';
options_.plot_shock_decomp.use_shock_groups = 'Joao';
options_.plot_shock_decomp.plot_end_date = dates('2020Q4');
options_.plot_shock_decomp.plot_init_date = dates('2018Q4');
options_.plot_shock_decomp.graph_format = char('fig','eps');
var_list_ = char('LYOBS_EA','LC_EA','LI_EA','LE_EA');
plot_shock_decomposition(M_, oo_, options_, var_list_);

%% revision: difference between last vintage and previous one 
oo3 = load('F:\GM\SummerF19\2019_03_05_SF19_NTREND_TRest_newOpts_BFL_noEndoFN_FWD_gmmu2_VM_assmpt0_2step_SF19_fix_noM_SD\gemc_results','oo_');
oo_.realtime_shock_decomposition = struct_operator(oo2_.jrc.assmpt.realtime_shock_decomposition,oo3.oo_.jrc.assmpt.realtime_shock_decomposition,'-');
oo_.conditional_shock_decomposition = struct_operator(oo2_.jrc.assmpt.conditional_shock_decomposition,oo3.oo_.jrc.assmpt.conditional_shock_decomposition,'-');
oo_.realtime_conditional_shock_decomposition = struct_operator(oo2_.jrc.assmpt.realtime_conditional_shock_decomposition,oo3.oo_.jrc.assmpt.realtime_conditional_shock_decomposition,'-');
oo_.realtime_forecast_shock_decomposition =struct_operator(oo2_.jrc.assmpt.realtime_forecast_shock_decomposition,oo3.oo_.jrc.assmpt.realtime_forecast_shock_decomposition,'-');
options_ = set_default_plot_shock_decomposition_options(options_);
options_.plot_shock_decomp.detail_plot = 1;
options_.plot_shock_decomp.interactive = 1;
options_.plot_shock_decomp.nodisplay = 1;
options_.plot_shock_decomp.realtime = 1;
options_.plot_shock_decomp.steadystate = 1;
options_.plot_shock_decomp.vintage = 0;
options_.plot_shock_decomp.write_xls = 1;
options_.plot_shock_decomp.fig_name = 'assmpt_revision';
options_.plot_shock_decomp.type = 'aoa';
options_.plot_shock_decomp.use_shock_groups = 'Joao';
options_.plot_shock_decomp.plot_end_date = dates('2020Q4');
options_.plot_shock_decomp.plot_init_date = dates('2018Q4');
options_.plot_shock_decomp.graph_format = char('fig','eps');
var_list_ = char('LYOBS_EA','LC_EA','LI_EA','LE_EA');
plot_shock_decomposition(M_, oo_, options_, var_list_);
