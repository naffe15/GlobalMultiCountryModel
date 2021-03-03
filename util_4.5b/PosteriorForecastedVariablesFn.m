function out=PosteriorForecastedVariablesFn(nfrcst, M_,options_,oo_,estim_params_, addsteady, exnam, exval, parname, parval)
% global bayestopt_
% 
% % Load and transform data.
% transformation = [];
% if options_.loglinear && ~options_.logdata
%     transformation = @log;
% end
% xls.sheet = options_.xls_sheet;
% xls.range = options_.xls_range;
% 
% if ~isfield(options_,'nobs')
%     options_.nobs = [];
% end
% dataset = initialize_dataset(options_.datafile,options_.varobs,options_.first_obs,options_.nobs,transformation,options_.prefilter,xls);
% Y = dataset.data;
% gend = dataset.info.ntobs;
% data_index = dataset.missing.aindex;
% missing_value = dataset.missing.state;
% bayestopt_.mean_varobs = dataset.descriptive.mean';

ex=zeros(M_.exo_nbr, nfrcst);
if nargin==8,
    for j=1:length(exnam),
        ix=strcmp(exnam{j},cellstr(M_.exo_names));
        ex(ix,:)=exval(:,j);
    end
end

DirectoryName = CheckPath('metropolis',M_.dname);
smooth_file_list = dir([DirectoryName,filesep,'*_smooth*.mat']);
inno_file_list = dir([DirectoryName,filesep,'*_inno*.mat']);
tmp_file_list = dir([DirectoryName,filesep,'*_param*.mat']);
jfile=0;
for j=1:length(tmp_file_list),
    if isempty(strfind(tmp_file_list(j).name,'irf')),
        jfile=jfile+1;
        param_file_list(jfile)=tmp_file_list(j);
    end
end
clear tmp_file_list jfile,
jsmooth_file=1;
jinno_file=1;
jsmooth=0;
jinno=0;
stock_inno = load([DirectoryName,filesep,inno_file_list(1).name],'stock');
load([DirectoryName,filesep,smooth_file_list(1).name],'stock')
stock_yfcast=zeros(M_.endo_nbr,nfrcst+1,size(stock,3));
B=0;
disp(['MH Posterior forecast ...']);
hh = dyn_waitbar(0,['Posterior forecast ...']);
for j=1:length(param_file_list),
    load([DirectoryName,filesep,M_.fname,'_param',int2str(j),'.mat'],'stock_ys')
    stock_params=load([DirectoryName,filesep,M_.fname,'_param',int2str(j),'.mat'],'stock');
    B=B+size(stock_params.stock,1);
    for i=1:size(stock_params.stock,1),
        jsmooth = jsmooth+1;
        if jsmooth>size(stock,3),
            jsmooth=1;
            stock=stock_yfcast;
            save([DirectoryName,filesep,M_.fname,'_fcast',int2str(jsmooth_file),'.mat'],'stock')
            jsmooth_file = jsmooth_file+1;
            load([DirectoryName,filesep,M_.fname,'_smooth',int2str(jsmooth_file),'.mat'],'stock')
            stock_yfcast=zeros(M_.endo_nbr,nfrcst+1,size(stock,3));
        end
        ahat = squeeze(stock(:,:,jsmooth));
        jinno = jinno+1;
        if jinno>size(stock_inno.stock,3),
            jinno=1;
            jinno_file = jinno_file+1;
            stock_inno = load([DirectoryName,filesep,M_.fname,'_inno',int2str(jinno_file),'.mat'],'stock');
        end
        ehat = squeeze(stock_inno.stock(:,:,jinno));
        M_ = set_all_parameters(stock_params.stock(i,:)',estim_params_,M_);
%         [T,R,SteadyState] = dynare_resolve(M_, options_, oo_);
%         ahat=ahat(oo_.dr.order_var,:)-repmat(SteadyState(oo_.dr.order_var),1,gend);
        ahat=ahat(oo_.dr.order_var,:)-repmat(stock_ys(i,oo_.dr.order_var)',1,size(ahat,2));
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
        dsteady = stock_ys(i,:)'-SteadyState;
        yfcst=zeros(M_.endo_nbr,nfrcst+1);
        yfcst(:,1)=ahat(:,end)+dsteady(oo_.dr.order_var);
        for jfcast=1:nfrcst,
            yfcst(:,jfcast+1)=T*yfcst(:,jfcast)+R*ex(:,jfcast);
        end
%         stock_yfcast(:,:,jsmooth)=yfcst;
        if options_.loglinear
            stock_yfcast(oo_.dr.order_var,:,jsmooth) = yfcst+ ...
                repmat(log(stock_ys(i,oo_.dr.order_var)'),1,nfrcst+1);
        else
            stock_yfcast(oo_.dr.order_var,:,jsmooth) = yfcst+ ...
                repmat(stock_ys(i,oo_.dr.order_var)',1,nfrcst+1);
        end
%         [alphahat,etahat,epsilonhat,alphatilde,SS,trend_coeff,aK] = ...
%             DsgeSmoother(stock_params.stock(i,:)',gend,Y,data_index,missing_value);
%         max(max(ahat-alphahat)),
        dyn_waitbar(i/size(stock_params.stock,1),hh,['Parameters block ',int2str(j),'/',int2str(length(param_file_list))]);
    end
end
dyn_waitbar_close(hh);
stock=stock_yfcast;
save([DirectoryName,filesep,M_.fname,'_fcast',int2str(jsmooth_file),'.mat'],'stock')
% pm3(M_.endo_nbr,nfrcst+1,jsmooth_file,1200,'Projections',...
%     '',char('E_GY','E_GBY','E_TBYN'),M_.endo_names_tex,M_.endo_names,...
%     char('E_GY','E_GBY','E_TBYN'),'Projection',DirectoryName,'_fcast');
pm3(M_.endo_nbr,nfrcst+1,jsmooth_file,B,'Projections',...
    '',M_.endo_names(1:M_.orig_endo_nbr, :),M_.endo_names_tex,M_.endo_names,...
    M_.endo_names(1:M_.orig_endo_nbr, :),'Projection',DirectoryName,'_fcast');



out=1;