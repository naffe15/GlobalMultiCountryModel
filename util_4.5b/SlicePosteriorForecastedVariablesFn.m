function out=SlicePosteriorForecastedVariablesFn(nfrcst, M_,options_,oo_,estim_params_, addsteady, exnam, exval, parname, parval)
global bayestopt_

% Load and transform data.
transformation = [];
if options_.loglinear && ~options_.logdata
    transformation = @log;
end
xls.sheet = options_.xls_sheet;
xls.range = options_.xls_range;

if ~isfield(options_,'nobs')
    options_.nobs = [];
end
dataset = initialize_dataset(options_.datafile,options_.varobs,options_.first_obs,options_.nobs,transformation,options_.prefilter,xls);
Y = dataset.data;
gend = dataset.info.ntobs;
data_index = dataset.missing.aindex;
missing_value = dataset.missing.state;
bayestopt_.mean_varobs = dataset.descriptive.mean';

ex=zeros(M_.exo_nbr, nfrcst);
if nargin==8,
    for j=1:length(exnam),
        ix=strcmp(exnam{j},cellstr(M_.exo_names));
        ex(ix,:)=exval(:,j);
    end
end

DirectoryName = CheckPath('SLICE',M_.dname);
load([DirectoryName,filesep,M_.fname,'_slice'],'XSIM');
istart=floor(0.2*size(XSIM,1))+1;
XSIM=XSIM(istart:end,:);
B=size(XSIM,1);
stock=XSIM;
save([DirectoryName,filesep,M_.fname,'_param'],'stock');
B=size(XSIM,1);
MAX_jsmooth=6;
MAX_jinno=50;
naK = length(options_.filter_step_ahead);
nahead = max(options_.filter_step_ahead);
stock_yfcast=zeros(M_.endo_nbr,nfrcst+1,MAX_jsmooth);
stock_smooth=zeros(M_.endo_nbr,gend,MAX_jsmooth);
stock_inno=zeros(M_.exo_nbr,gend,MAX_jinno);
if naK,
    stock_filter_step_ahead=zeros(naK,M_.endo_nbr,gend+nahead,MAX_jsmooth);
end
jsmooth_file = 1;
jinno_file = 1;
jsmooth = 0;
jinno = 0;
disp(['SLICE Posterior smoother and forecast ...']);
hh = dyn_waitbar(0,['Slice smoother and forecast ...']);
for i=1:size(XSIM,1),
    jsmooth = jsmooth+1;
    if jsmooth>MAX_jsmooth,
        jsmooth=1;
        if nfrcst,
        stock=stock_yfcast;
        save([DirectoryName,filesep,M_.fname,'_fcast',int2str(jsmooth_file)],'stock')
        stock_yfcast=zeros(M_.endo_nbr,nfrcst+1,MAX_jsmooth);
        end
        stock=stock_smooth;
        save([DirectoryName,filesep,M_.fname,'_smooth',int2str(jsmooth_file)],'stock')
        stock_smooth=zeros(M_.endo_nbr,gend,MAX_jsmooth);
        if naK
            stock=stock_filter_step_ahead;
            save([DirectoryName,filesep,M_.fname,'_filter_step_ahead',int2str(jsmooth_file)],'stock')
            stock_filter_step_ahead=zeros(naK,M_.endo_nbr,gend+nahead,MAX_jsmooth);
        end
        jsmooth_file = jsmooth_file+1;
    end
    jinno = jinno+1;
    if jinno>MAX_jinno,
        jinno=1;
        stock=stock_inno;
        save([DirectoryName,filesep,M_.fname,'_inno',int2str(jinno_file)],'stock')
        stock_inno=zeros(M_.exo_nbr,gend,MAX_jinno);
        jinno_file = jinno_file+1;
    end
    M_ = set_all_parameters(XSIM(i,:)',estim_params_,M_);
    [T,R,SteadyState] = dynare_resolve(M_, options_, oo_);
    [alphahat,etahat,epsilonhat,alphatilde,SS,trend_coeff,aK] = ...
        DsgeSmoother(XSIM(i,:)',gend,Y,data_index,missing_value);
    if nfrcst,
    if nargin==10,
        for jpar=1:length(parname),
            ip = strmatch(parname{jpar},M_.param_names,'exact');
            if isempty(ip)
                error(['Parameter name ' parname{jpar} ' doesn''t exist'])
            end
            M_.params(ip) = parval(jpar);
        end
    end
    [T,R,SteadyState] = dynare_resolve(M_, options_, oo_);
    dsteady = SS-SteadyState;
    yfcst=zeros(M_.endo_nbr,nfrcst+1);
    yfcst(:,1)=alphahat(:,end)+dsteady(oo_.dr.order_var);
    for jfcast=1:nfrcst,
        yfcst(:,jfcast+1)=T*yfcst(:,jfcast)+R*ex(:,jfcast);
    end
    end
    %         stock_yfcast(:,:,jsmooth)=yfcst;
    stock_inno(:,:,jinno) = etahat;
    if naK
        stock_filter_step_ahead(:,oo_.dr.order_var,:,jsmooth) = aK(options_.filter_step_ahead,1:M_.endo_nbr,:);
    end
    if options_.loglinear
        if nfrcst,
        stock_yfcast(oo_.dr.order_var,:,jsmooth) = yfcst+ ...
            repmat(log(SS(oo_.dr.order_var)),1,nfrcst+1);
        end
        stock_smooth(oo_.dr.order_var,:,jsmooth) = alphahat+ ...
            repmat(log(SS(oo_.dr.order_var)),1,gend);
    else
        if nfrcst,
        stock_yfcast(oo_.dr.order_var,:,jsmooth) = yfcst+ ...
            repmat(SS(oo_.dr.order_var),1,nfrcst+1);
        end
        stock_smooth(oo_.dr.order_var,:,jsmooth) = alphahat+ ...
            repmat(SS(oo_.dr.order_var),1,gend);
    end
    %         max(max(ahat-alphahat)),
    dyn_waitbar(i/B,hh);
end
dyn_waitbar_close(hh);
if nfrcst,
stock=stock_yfcast(:,:,1:jsmooth);
save([DirectoryName,filesep,M_.fname,'_fcast',int2str(jsmooth_file)],'stock')
pm3(M_.endo_nbr,nfrcst+1,jsmooth_file,B,'Projections',...
    '',M_.endo_names(1:M_.orig_endo_nbr, :),M_.endo_names_tex,M_.endo_names,...
    M_.endo_names(1:M_.orig_endo_nbr, :),'Projection',DirectoryName,'_fcast');
end
stock=stock_smooth(:,:,1:jsmooth);
save([DirectoryName,filesep,M_.fname,'_smooth',int2str(jsmooth_file)],'stock')
stock=stock_inno(:,:,1:jinno);
save([DirectoryName,filesep,M_.fname,'_inno',int2str(jinno_file)],'stock')
if naK
stock=stock_filter_step_ahead(:,:,:,1:jsmooth);
save([DirectoryName,filesep,M_.fname,'_filter_step_ahead',int2str(jsmooth_file)],'stock')
end
pm3(M_.endo_nbr,gend,jsmooth_file,B,'Smoothed variables',...
    '',M_.endo_names(1:M_.orig_endo_nbr, :),M_.endo_names_tex,M_.endo_names,...
    M_.endo_names(1:M_.orig_endo_nbr, :),'SmoothedVariables',DirectoryName,'_smooth');
pm3(M_.exo_nbr,gend,jinno_file,B,'Smoothed shocks',...
    '',M_.exo_names,M_.exo_names_tex,M_.exo_names,...
    M_.exo_names,'SmoothedShocks',DirectoryName,'_inno');


out=1;