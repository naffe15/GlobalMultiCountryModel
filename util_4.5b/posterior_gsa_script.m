%% load results
clear all
clc
dynare_root=dynare_config();
load gemc_results
load data T

%% get the samples
varlist = char('GYOBS_EA','PHIYOBS_EA','INOM_EA','GE_EA');
%varlist = char('YGAP_EA','NTREND_EA');
varlist = [];
[ oo, params, varlist ] = get_posterior_sample( 'forecast', varlist, M_ ); % loads posterior draws of forecasts in oo and associated parameters in params
[ oo1, params1, varlist ] = get_posterior_sample( 'smooth', varlist, M_ ); % loads posterior draws of smoother in oo1 and associated parameters in params1
[ oo2, params2, varlist ] = get_posterior_sample( 'irf', varlist, M_ ); % loads posterior draws of irfs in oo2 and associated parameters in params2
[ oo3, params3, varlist ] = get_posterior_sample( 'filter', varlist, M_ ); % loads posterior draws of smoofilter steps ahead  in oo3 structure and associated parameters in params1
[ lpost, params0, parlist ] = get_posterior_sample('param', [], M_, bayestopt_ ); % loads posterior params in params0


%% do irf plots
figure,
for j=1:size(oo2.e_A,1)
    subplot(3,4,j)
    plot(squeeze(oo2.e_A(j,:,:))),
    title(varlist(j,:),'interpreter','none')
end

%% visual irf MC
ind_t = 6;
figure, scatter_plots(squeeze(oo2.e_A(1:3,ind_t,:))',params2, varlist(1:3,:),[], [], [], [], [], options_);


%% do plots forecasts
figure,
for j=1:size(oo,1)
    subplot(2,2,j)
    plot(T(options_.first_obs:end),cat(1,squeeze(oo1(j,:,:)),squeeze(oo(j,:,:)))),
    title(varlist(j,:),'interpreter','none')
end

%% visual MC
ind_t = 6;
figure, scatter_plots(squeeze(oo(:,ind_t,:))',params, varlist,[], [], [], [], [], options_);

%% select run
figure,
irun = 57;
for j=1:size(oo,1)
    subplot(2,2,j)
    plot(T(options_.first_obs:end),cat(2,squeeze(oo1(j,:,irun)),squeeze(oo(j,:,irun)))),
    title(varlist(j,:),'interpreter','none')
end

%% do plots smoother
figure,
for j=1:size(oo1,1)
    subplot(2,2,j)
    plot(T(options_.first_obs:end-options_.forecast),squeeze(oo1(j,:,:))),
    title(varlist(j,:),'interpreter','none')
end


%% MCF setting
addpath z:\svn\mcf\matlab\
options_mcf.pvalue_corr=0;
options_mcf.param_names=bayestopt_.name;


%% filter lowest 20% for inom
nsam = size(params,1);
[s, is] = sort(oo(3,ind_t,:));
ibeha = is(1:ceil(nsam*0.2));
ino = is(ceil(nsam*0.2)+1:end);
inx=mcf_analysis(params, ibeha, ino,options_mcf);

figure,
for j=1:size(oo,1)
    subplot(2,2,j)
    plot(T(options_.first_obs:end),cat(1,squeeze(oo1(j,:,ibeha)),squeeze(oo(j,:,ibeha)))),
    title(varlist(j,:),'interpreter','none')
end

%% filter highest 20% for GE
ind_t = 1;
nsam = size(params,1);
[s, is] = sort(-oo(4,ind_t,:));
ibeha = is(1:ceil(nsam*0.2));
ino = is(ceil(nsam*0.2)+1:end);
inx=mcf_analysis(params, ibeha, ino,options_mcf);

[s, is] = sort(-oo(3,ind_t,:));

%% plot posterior densities
load gemc_results_slice mode1
load(mode1,'xparam1')
for j=1:length(inx),
    
    pnam = bayestopt_.name{inx(j)};
    subplot(3,3,j)
    try
        plot(oo_.posterior_density.parameters.(pnam)(:,1),oo_.posterior_density.parameters.(pnam)(:,2),'linewidth',2)
    catch
        plot(oo_.posterior_density.shocks_std.(pnam)(:,1),oo_.posterior_density.shocks_std.(pnam)(:,2),'linewidth',2)
    end
    yl=get(gca,'ylim');
    hold all, plot(xparam1([inx(j) inx(j)]),yl,'linewidth',2,'linestyle','--')
    set(gca,'ylim',yl);
    title(pnam,'interpreter','none')
    
end

%% compare param estimates
parlist = cellstr(char('SIGMAZ_EA','SIGMAZ_US','SIGMAC_RoW'));
parlist = cellstr(char('RHO_PX_US_RoW','RHO_PX_EA_RoW'));
parlist = cellstr(char('RHO_BW_RoW','RHO_BW_EA'));
figure, 
for j=1:size(parlist,1),
    try
        plot(oo_.posterior_density.parameters.(parlist{j})(:,1),oo_.posterior_density.parameters.(parlist{j})(:,2),'linewidth',2)
    catch
        plot(oo_.posterior_density.shocks_std.(parlist{j})(:,1),oo_.posterior_density.shocks_std.(parlist{j})(:,2),'linewidth',2)
    end
    hold all,
end
hl=legend(parlist);
set(hl,'interpreter','none')

%%
addpath z:\svn\mcf\matlab\
%%
parlist = cellstr(char('RHO_BW_RoW','RHO_BW_EA'));
parlist = cellstr(char('ETAIPHI_EA','RHO_BW_EA'));
parlist = cellstr(char('ETAIPHI_RoW','RHO_BW_RoW'));
parlist = cellstr(char('ETAIPHI_EA','ETAIPHI_RoW', 'ETAIPHI_US'));
parlist = cellstr(char('RHOINOM_EA','RHOINOM_RoW', 'RHOINOM_US'));
parlist = cellstr(char('RHOINOM_EA','RHO_BW_EA'));
parlist = cellstr(char('RHOINOM_RoW','RHO_BW_RoW'));
parlist = cellstr(char('GAMMAWR_EA','RHO_TZEPS_MUY_EA','T_MUY_Y0_EA'));

parlist = cellstr(char('ALPHA_SLI_EA','AS1_EA','SLI_SS_EA','GAMMAI1_EA','GAMMAI2_EA'));
pindx = find(ismember(bayestopt_.name,parlist));
pscatter(params0(:,pindx))
%%
figure, 
subplot(221)
plot(reshape(params0(:,pindx(1)),50,4)), 
title(bayestopt_.name(pindx(1)),'interpreter','none')
subplot(222)
plot(reshape(params0(:,pindx(2)),50,4)), 
title(bayestopt_.name(pindx(2)),'interpreter','none')
subplot(223)
plot(reshape(params0(:,pindx(3)),50,4)), 
title(bayestopt_.name(pindx(3)),'interpreter','none')
% subplot(224)
% plot(reshape(params0(:,pindx(4)),50,4)), 
% title(bayestopt_.name(pindx(4)),'interpreter','none')


%% FIT
pplotvar = 'GYOBS_EA';
TT = T(options_.first_obs:options_.first_obs+options_.nobs-1);
figure
    plot(TT,oo_.UpdatedVariables.(pplotvar),'linewidth',2)
    hold on,
    forecast = size(oo3.(pplotvar),1);
    np = size(oo3.(pplotvar),3);
    for jt=1:length(oo_.UpdatedVariables.(pplotvar))
        for jp=1:np
            qforecastedvar(jp,1) = oo_.UpdatedVariables.(pplotvar)(jt);
            qforecastedvar(jp,2:forecast+1)=diag(squeeze(oo3.(pplotvar)(:,jt+1:jt+forecast,jp)));
        end
        plot(TT(jt):(TT(jt)+forecast),qforecastedvar)
    end
    plot(TT,oo_.UpdatedVariables.(pplotvar),'linewidth',2,'color',[0 0 0])
    hold off,
