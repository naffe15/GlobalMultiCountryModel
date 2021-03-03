function [moments, data_moments] = check_covariances(do_the_plot,var_plot,varargin)

% function check_covariances(varargin);
% compares model theoretical  moments with those in the data
% (using SmoothedVariables)
% ACF and ACF for the data are computed without  prefiltering, i.e.
% consitent with model estimation where data fluctuate around the  model
% steady state and not the sample mean
%
% inputs: list of endogenous
%
% outputs:
% moments: wrap-up of model implied ACF/CCF for the declared list of
%     edogenous
% data_moments: wrap-up of data ACF/CCF for the declared list of
%     edogenous
%
global M_ oo_ options_

first_obs = options_.presample+1;

lambda=options_.hp_filter;
% [s,desvabs] = hpfilter(y,w,plotter);

fname_ = M_.fname;
[DirectoryName, info] = CheckPath('Output',fname_);
% //oo1 = oo_;
% //moments.gamma_y = oo_.gamma_y;
moments.mean = oo_.mean;
moments.var = oo_.var;
moments.corr = oo_.var./(sqrt(diag(oo_.var))*sqrt(diag(oo_.var))');
nlags = options_.ar;
for j=1:nlags
    moments.autocorr(:,j) = diag(oo_.autocorr{j});
end

moments.acf=struct();
moments.ccf=struct();
data_moments.var=zeros(length(varargin),length(varargin));
data_moments.corr=zeros(length(varargin),length(varargin));
data_moments.autocorr=zeros(length(varargin),nlags);
data_moments.acf=struct();
data_moments.ccf=struct();
for j=1:length(varargin)
    moments.acf.(varargin{j}) = moments.autocorr(j,:)';
end
if ~isempty(options_.datafile)
    SmoothedVariables = load(options_.datafile);
    fnames = fieldnames(SmoothedVariables);
    for k=1:length(fnames)
        SmoothedVariables.(fnames{k}) = SmoothedVariables.(fnames{k})(options_.first_obs:options_.first_obs+options_.nobs-1);
    end
end
SmoothedVariables0 = [];
if isfield(oo_, 'SmoothedVariables')
    SmoothedVariables0 = oo_.SmoothedVariables;
end
nobs = length(SmoothedVariables.(fnames{1}));
nsiz = length(SmoothedVariables.(fnames{1}))-first_obs+1;
nvar = length(varargin);
data = NaN(nsiz, nvar);
cc = NaN(nvar,nvar);
cova = NaN(nvar,nvar);
std_dev = NaN(nvar,1);
for j=1:nvar
    if isfield(SmoothedVariables,varargin{j})
        tmp = SmoothedVariables.(varargin{j})-get_mean(varargin{j});
    elseif ~isempty(SmoothedVariables0)
        tmp = SmoothedVariables0.(varargin{j})-get_mean(varargin{j});
    else        
        tmp = NaN(nobs,1);
    end
    data(:,j) = tmp(first_obs:end);
    if lambda
        [s,desvabs] = hpfilter(data(:,j),lambda);
        data(:,j) = data(:,j) - s;
    end
    cova(j,j)=sum(data(:,j).^2)/size(data,1);
    std_dev(j) = sqrt(cova(j,j));
    cc(j,j)=1;
    for ii=1:j-1
        cova(j,ii)=sum(data(:,j).*data(:,ii))/size(data,1);
        cova(ii,j)=cova(j,ii);
        cc(j,ii)=cova(j,ii)/std_dev(j)/std_dev(ii);
        cc(ii,j) = cc(j,ii);
    end
end

if lambda
    hp_filt_txt='_hpfilt';
else
    hp_filt_txt='';
end
data_moments.var=cova;
data_moments.corr=cc;
n=size(data,1);
if do_the_plot
    nvar=size(data,2);
    nrow=ceil(sqrt(nvar));
    if nrow*(nrow-1)>=nvar
        ncol=nrow-1;
    else
        ncol=nrow;
    end
    hACF = dyn_figure(options_.nodisplay,'Name',['ACF - model vs. data']);
    posit = get(0,'monitorposition');
    posit=posit(end,:);
end
for ii=-nlags:nlags
    sigmaXCF(ii+nlags+1,1) = sqrt(1/(n-abs(ii)));
    BoundsXCF(ii+nlags+1,[1 2]) = sigmaXCF(ii+nlags+1,1)*[2 -2];
end
count = -1;
%icount = 0;
for j=1:length(varargin)
    for j2=1:length(var_plot)
        indj2 = strmatch(var_plot{j2},varargin, 'exact');
        if indj2==j
            count=count+1;
        end
    end
    % autocorr/crosscorr: 0: prefilter: assume data are without the mean,
    % consistent with moments computed around the steady state of  the
    % model
    [ACF,Lags] = autocorr(data(:,j),nlags,0,2,0);
    for ii=1:nlags
        sigmaQ(ii,1) = sqrt((1 + 2*(ACF(2:ii)'*ACF(2:ii)))/n);
        Bounds(ii,[1 2]) = sigmaQ(ii,1)*[2 -2];
    end
    data_moments.autocorr(j,:)=ACF(2:end);
    data_moments.acf.(varargin{j}) = [Lags(2:end) ACF(2:end) Bounds];
    %
    %  Plot the sample ACF.
    %
    %    figure(hACF)
    if do_the_plot
        set(0,'CurrentFigure',hACF)
        subplot(nrow,ncol,j)
        lineHandles  =  stem([1:nlags]-0.1,moments.autocorr(j,:), 'filled' , 'k-o');
        set   (lineHandles(1) , 'MarkerSize' , 2, 'MarkerFaceColor', 'k')
        xlabel('Lag','fontsize',8)
        ylabel('acf','fontsize',8)
        title (varargin{j},'interpreter','none','fontsize',8)
        hold  ('on')
        %
        %  Plot the confidence bounds under the hypothesis that the underlying
        %  Series is really an MA(Q) process. Bartlett's approximation gives
        %  an indication of whether the ACF is effectively zero beyond lag Q.
        %  For this reason, the confidence bounds (horizontal lines) appear
        %  over the ACF ONLY for lags GREATER than Q (i.e., Q+1, Q+2, ... nLags).
        %  In other words, the confidence bounds enclose ONLY those lags for
        %  which the null hypothesis is assumed to hold.
        %
        
        lineHandles  =  stem([1:nlags] , ACF(2:end),'filled' , 'b-o');
        set   (lineHandles(1) , 'MarkerSize' , 2)
        plot([1:nlags] , Bounds , '--r');
        plot([1 nlags] , [0 0] , '-b');
        hold('off')
        a  =  axis;
        axis([0.5 nlags a(3) 1]);
    end
    for j2=1:length(var_plot)
        indj2 = strmatch(var_plot{j2},varargin, 'exact');
        if indj2==j && do_the_plot
            hCCF = dyn_figure(options_.nodisplay,'Name',[varargin{j} ' - CCF''s - model vs. data']);
        end
    end
    
    for i=1:length(varargin)
        %       if count == 0
        %           icount = icount+1;
        %       end
        [XCF,Lags] = crosscorr(data(:,j),data(:,i),nlags,2,0);
        
        for k=-nlags:nlags
            if k<0
                tmp1(k+nlags+1,1)=oo_.autocorr{abs(k)}(i,j);
            elseif  k==0
                tmp1(k++nlags+1,1)=moments.corr(i,j);
            else
                tmp1(k++nlags+1,1)=[oo_.autocorr{abs(k)}(j,i)];
            end
        end
        
        
        
        %
        %  Plot the sample XCF.
        %
        
        if do_the_plot
            if isempty(var_plot)
                var_plot = varargin;
            end
            for j2=1:length(var_plot)
                indj2 = strmatch(var_plot{j2},varargin, 'exact');
                if indj2==j
                    %lastst = var_plot{j2}(end-2:end);
                    %indi2 = strmatch(lastst,varargin{i}(end-2:end),'exact');
                    %if indi2 ==1,
                    %nsub = i-count*(icount);
                    %subplot(nrow,ncol,nsub)
                    subplot(nrow,ncol,i)
                    lineHandles  =  stem([-nlags:nlags]-0.1,tmp1, 'filled' , 'k-o');
                    set   (lineHandles(1) , 'MarkerSize' , 2, 'MarkerFaceColor', 'k')
                    xlabel('lag','fontsize',8)
                    ylabel('ccf','fontsize',8)
                    %         text (-nlags-0.5,1.1,[varargin{j} '(t)         ' varargin{i} '(t-lag)'],'interpreter','none','fontsize',8)
                    %         text (0.1,1.1,[varargin{i} '(t-lag)'],'interpreter','none','fontsize',8,'units','normalized')
                    title([varargin{i} '(t-lag)'],'interpreter','none','fontsize',8)
                    hold  ('on')
                    %
                    %  Plot the confidence bounds under the hypothesis
                    %  that the underlying series are uncorrelated.
                    %
                    lineHandles  =  stem(Lags , XCF , 'filled' , 'b-o');
                    set   (lineHandles(1) , 'MarkerSize' , 2)
                    plot([-nlags:nlags] , BoundsXCF , '--r');
                    
                    plot([-nlags nlags] , [0 0] , '-b');
                    hold('off')
                    a  =  axis;
                    axis([-nlags-0.5 nlags+0.5 a(3:end)]);
                    %end
                end
            end
            
        end
        data_moments.ccf.([varargin{j} '_' varargin{i}]) = [Lags XCF BoundsXCF];
        moments.ccf.([varargin{j} '_' varargin{i}]) = [Lags tmp1];
    end
    for j2=1:length(var_plot)
        indj2 = strmatch(var_plot{j2},varargin, 'exact');
        if indj2==j
            if do_the_plot
                
                if ~isoctave
                    annotation('textbox', [0.1,0.02,0.35,0.05],'String', 'DATA','Color','Blue','horizontalalignment','center','interpreter','none');
                    annotation('textbox', [0.55,0.02,0.35,0.05],'String', 'MODEL','Color','Black','horizontalalignment','center','interpreter','none');
                end
                set(gcf, ...
                    'position',posit, ...
                    'paperorientation','landscape', ...
                    'papertype','a3', ...
                    'PaperPositionMode','auto')
                dyn_saveas(hCCF,[DirectoryName,filesep,fname_,'_CCF',hp_filt_txt,'_',varargin{j}],options_.nodisplay,options_.graph_format)
            end
        end
    end
end
if do_the_plot
    
    if ~isoctave
        annotation('textbox', [0.1,0.02,0.35,0.05],'String', 'DATA','Color','Blue','horizontalalignment','center','interpreter','none');
        annotation('textbox', [0.55,0.02,0.35,0.05],'String', 'MODEL','Color','Black','horizontalalignment','center','interpreter','none');
    end
    set(hACF, ...
        'position',posit, ...
        'paperorientation','landscape', ...
        'papertype','a3', ...
        'PaperPositionMode','auto')
    dyn_saveas(hACF,[DirectoryName,filesep,fname_,'_ACF',hp_filt_txt],options_.nodisplay,options_.graph_format)
end
end




% data_moments.var =  out1;
% data_moments.var =  cov_nan([E_Y E_C E_I E_WRPHI E_L E_PHIC E_G]);
%
% data_moments.corr =  corrcoef_nan([E_GY E_GC E_GI E_WRPHI E_GL E_PHI E_PHIC E_GG]);
% data43 = [E_GY E_GC E_GI E_WRPHI E_GL E_PHI E_PHIC E_GG ];
% for j=1:size(data43,2),
% datx=data43(:,j);
% [acov(j,:), acorr(j,:)] = sample_autocovariance(datx(~isnan(datx)),5);
% end
% data_moments.autocorr = acorr(:,2:end);
