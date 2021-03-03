function [XSIM,dataset_,xparam1, M_, options_, oo_, estim_params_,bayestopt_] = dyn_slice_loop(G, M_,options_,oo_,estim_params_,xstart)

if ~isfield(M_,'dname')
    M_.dname=M_.fname;
end
MhDirectoryName = CheckPath('SLICE',M_.dname);
ModelName = M_.fname;
drop=options_.mh_drop;

[dataset_, dataset_info, xparam1, hh, M_, options_, oo_, estim_params_, bayestopt_, bounds] = ...
    dynare_estimation_init([], M_.dname, 0, M_, options_, oo_, estim_params_);
iload=0;
if nargin ==6,
    if strcmp(xstart,'load_slice_file'),
        load([MhDirectoryName '/' ModelName '_slice.mat'],'XSIM','logpo2','NEVAL');
        iload=1;
        xparam1=XSIM(end,:)';
        XSIM0=XSIM;
        NEVAL0=NEVAL;
        logpo20=logpo2;
        G0=size(XSIM,1);
        save([MhDirectoryName '/' ModelName '_slice0.mat'],'XSIM','logpo2','NEVAL');
    else
        xparam1=xstart;
    end
end
XSIM = zeros(G,2);
start= tic;
thetaprior = [xparam1 xparam1 bounds.lb bounds.ub];
nt = length(xparam1);
theta = xparam1;
h = dyn_waitbar(0,'SLICE loop ...');
tot_NEVAL=0;
[record.InitialSeeds.Unifor,record.InitialSeeds.Normal] = get_dynare_random_generator_state();
save([MhDirectoryName '/' ModelName '_slice_history.mat'],'theta','record');
for g = 1:G    
    dyn_waitbar((g-1)/G,h,['Performing SLICE iteration ', int2str(g),'/',int2str(G)]);
    for it = 1:nt
        [NEVAL,XSIM(g,it)] = slice_dynare(it,theta,thetaprior(it,:),dataset_,dataset_info, options_, M_, estim_params_,bayestopt_,bounds, oo_);
        theta(it) = XSIM(g,it);
        tot_NEVAL=tot_NEVAL + NEVAL;
    end
    logpo2(g)     = -dsge_likelihood(theta,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_);
    save([MhDirectoryName '/' ModelName '_slice.mat'],'XSIM','logpo2','NEVAL');
end
[record.LastSeeds.Unifor,record.LastSeeds.Normal] = get_dynare_random_generator_state();
save([MhDirectoryName '/' ModelName '_slice_history.mat'],'record','-append');
dyn_waitbar_close(h);
if iload, 
    XSIM = [XSIM0;XSIM];
    logpo2 = [logpo20 logpo2];
    NEVAL = NEVAL+NEVAL0;
    G=size(XSIM,1);
    copyfile([MhDirectoryName '/' ModelName '_slice.mat'],[MhDirectoryName '/' ModelName '_slice1.mat']);
    save([MhDirectoryName '/' ModelName '_slice.mat'],'XSIM','logpo2','NEVAL');
end
xparam1=mean(XSIM(drop*G+1:end,:))';
SIGMA=cov(XSIM(drop*G+1:end,:));
hh=inv(SIGMA);
parameter_names = bayestopt_.name;
save([M_.fname,'_mean.mat'],'xparam1','hh','SIGMA','parameter_names')
[ml,im]=max(logpo2);
xparam1=XSIM(im,:)';
save([M_.fname,'_mode_slice.mat'],'xparam1','hh','SIGMA','parameter_names')
if exist([M_.fname,'_mode.mat'],'file'),
    tmp = load([M_.fname,'_mode.mat'],'parameter_names');
    if isequal(tmp.parameter_names,parameter_names),
        load([M_.fname,'_mode.mat'],'xparam1')
    else
        skipline,
        disp('Warning')
        disp(['The mode file ' M_.fname,'_mode.mat has different prior list w.r.t. current SLICE loop!'])
        skipline,
    end
end
nvx     = estim_params_.nvx;
nvn     = estim_params_.nvn;
ncx     = estim_params_.ncx;
ncn     = estim_params_.ncn;
np      = estim_params_.np ;
DirectoryName = CheckPath('Output',M_.dname);

xparam=xparam1;
ifig=0;
for indx=1:nt,
    if mod(indx,9)==1,
        ifig=ifig+1;
        h=dyn_figure(options_.nodisplay,'Name','SLICE Priors and Posteriors');
        iplo=0;
    end
    iplo=iplo+1;
    [post_mean, post_median, post_var, hpd_interval, post_deciles, ...
        density] = posterior_moments(XSIM(drop*G+1:end,indx),1,options_.mh_conf_sig);
    name = bayestopt_.name{indx};
    if indx>(nvx+nvn+ncx+ncn),
        type = 'parameters';
    elseif indx<= nvx,
        type = 'shocks_std';
    elseif indx<= nvx+nvn,
        type = 'measurement_errors_std';
    elseif indx<= nvx+nvn+ncx,
        type = 'shocks_corr';
    elseif indx<= nvx+nvn+ncx+ncn,
        type = 'measurement_errors_corr';
    end
    oo_ = Filloo(oo_,name,type,post_mean,hpd_interval,post_median,post_var,post_deciles,density);
    [x,f,abscissa,dens,binf,bsup] = draw_prior_density(indx,bayestopt_);
    subplot(3,3,iplo),
    [am,im]=max(density(:,2));
    xparam(indx,1)=density(im,1);
    plot(abscissa,dens,'LineWidth',2,'Color',[0.7 0.7 0.7])
    hold all, plot(density(:,1),density(:,2), 'b','LineWidth',2)
    yl=get(gca,'ylim');
    set(gca,'ylim', [0 yl(2)]),
    set(gca,'xlim', [min(abscissa(1),density(1,1)) max(abscissa(end),density(end,1))]),
%     set(gca,'xlim', [bayestopt_.lb(indx) bayestopt_.ub(indx)]),
    plot(xparam1([indx indx]),[0 yl(2)],'g--','LineWidth',2)
    title(bayestopt_.name{indx},'interpreter','none')
    if iplo==9 || indx==nt,
        dyn_saveas(h,[DirectoryName '/' M_.fname '_PriorsAndPosteriors_SLICE' int2str(ifig)],options_.nodisplay,options_.graph_format);
    end
end
xparam1=xparam;
llo  = -dsge_likelihood(xparam,dataset_, dataset_info, options_,M_,estim_params_,bayestopt_,bounds,oo_);
if llo>ml
    save([M_.fname,'_mode_slice.mat'],'xparam1','-append')
end
save([M_.fname,'_mean.mat'],'xparam','-append')
return

function oo = Filloo(oo,name,type,postmean,hpdinterval,postmedian,postvar,postdecile,density)
eval(['oo.posterior_mean.' type '.' name ' = postmean;']);
eval(['oo.posterior_hpdinf.' type '.' name ' = hpdinterval(1);']); 
eval(['oo.posterior_hpdsup.' type '.' name ' = hpdinterval(2);']);      
eval(['oo.posterior_median.' type '.' name ' = postmedian;']);
eval(['oo.posterior_variance.' type '.' name ' = postvar;']);
eval(['oo.posterior_deciles.' type '.' name ' = postdecile;']);
eval(['oo.posterior_density.' type '.' name ' = density;']);

function [post_mean,hpd_interval,post_var] = Extractoo(oo,name,type)
hpd_interval = zeros(2,1);
eval(['post_mean = oo.posterior_mean.' type '.' name ';']);
eval(['hpd_interval(1) = oo.posterior_hpdinf.' type '.' name ';']); 
eval(['hpd_interval(2) = oo.posterior_hpdsup.' type '.' name ';']);
eval(['post_var = oo.posterior_variance.' type '.' name ';']);