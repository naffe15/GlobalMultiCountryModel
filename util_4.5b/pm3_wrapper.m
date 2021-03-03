pm3(M_.endo_nbr,dataset_.nobs,4,options_.B,'Smoothed variables',...
        '','CU_IT',M_.endo_names_tex,M_.endo_names,...
'CU_IT','SmoothedVariables',CheckPath('metropolis',M_.dname),'_smooth'); 


%% this is for irfs

plot_IRF_post(xname2, [var_irf], oo_.PosteriorIRF.dsge,[],4,3,'Output_test')


%% this is for shock decompositions
posterior_shock_decomp_smooth_q(M_, options_, oo_, estim_params_,2000, ex_names_, leg, vname)
posterior_shock_decomp_smooth_2(M_, options_, oo_, estim_params_,2000, ex_names_, leg, {'GYOBS_IT','PHIYOBS_IT'})