function out=CrossCovariance(M_, oo_, options_, iplot)
% out=CrossCovariance(M_, oo_, options_, iplot)
% computes moments, ACF CCF of smoothed shocks, and optionally plots them

[DirectoryName, info] = CheckPath('Output',M_.fname);

%if nargin>2
%    first_obs = options_.presample+1;
%else
    first_obs = 1;
%end

if nargin<4
    iplot=0;
end
% check covariance of shocks
ccova = NaN(M_.exo_nbr,M_.exo_nbr);
cc = NaN(M_.exo_nbr,M_.exo_nbr);
std_dev = NaN(M_.exo_nbr,1);
eps1 = NaN(length(oo_.SmoothedShocks.(deblank(M_.exo_names(1,:))))-first_obs+1,M_.exo_nbr);
for j=1:M_.exo_nbr
    tmp = oo_.SmoothedShocks.(deblank(M_.exo_names(j,:)));
    eps1(:,j)= tmp(first_obs:end);
    ccova(j,j)=sum(eps1(:,j).^2)/size(eps1,1);
    %     ccova(j,j)=sum((eps1(:,j)-mean(eps1(:,j))).^2)/size(eps1,1);
    std_dev(j) = sqrt(ccova(j,j));
    cc(j,j)=1;
    for ii=1:j-1
        %         ccova(j,ii)=sum((eps1(:,j)-mean(eps1(:,j))).*(eps1(:,ii)-mean(eps1(:,ii))))/size(eps1,1);
        ccova(j,ii)=sum((eps1(:,j)).*(eps1(:,ii)))/size(eps1,1);
        ccova(ii,j)=ccova(j,ii);
        if std_dev(j) && std_dev(ii)
            cc(j,ii)=ccova(j,ii)/std_dev(j)/std_dev(ii);
            cc(ii,j)=cc(j,ii);
        end
    end
end

% cc = corrcoef(eps1);
% ccova = cov(eps1);
out.cc=cc;
out.ccova=ccova;

xls_cell=num2cell(ccova);
xls_cell=[cellstr(M_.exo_names) xls_cell];
xls_cell=[cellstr(char('',M_.exo_names))'; xls_cell];
if ~ismac
    [SUCCESS,MESSAGE]=xlswrite([M_.fname '_moments_shocks.xls'],xls_cell,'covariance');
else
%     [SUCCESS,MESSAGE]=xlswrite_MACOS([M_.fname '_moments_shocks.xls'],xls_cell,'covariance');
        writetable(array2table(xls_cell), [M_.fname '_moments_shocks.xls'], 'Sheet', 'covariance'); 
end

xls_cell=num2cell(cc);
xls_cell=[cellstr(M_.exo_names) xls_cell];
xls_cell=[cellstr(char('',M_.exo_names))'; xls_cell];
if ~ismac
    [SUCCESS,MESSAGE]=xlswrite([M_.fname '_moments_shocks.xls'],xls_cell,'correlation');
else
%     [SUCCESS,MESSAGE]=xlswrite_MACOS([M_.fname '_moments_shocks.xls'],xls_cell,'correlation');
        writetable(array2table(xls_cell), [M_.fname '_moments_shocks.xls'], 'Sheet', 'correlation'); 
end
posit = get(0,'monitorposition');
posit=posit(end,:);

names=cellstr(M_.exo_names);
variables_=names;
if  iplot
    fname_ = M_.fname;
    nlags=5;
    indx1=find(sqrt(diag(ccova)));
    data=eps1(:,indx1);
    variables_=cellstr(M_.exo_names(indx1,:));
    n=size(data,1);
    nvar=size(data,2);
    nrow=ceil(sqrt(nvar));
    if nrow*(nrow-1)>=nvar,
        ncol=nrow-1;
    else
        ncol=nrow;
    end
    
    out.autocorr=zeros(length(variables_),nlags);
    out.acf=struct();
    out.ccf=struct();
    
    hACF = dyn_figure(options_.nodisplay,'Name',['ACF - shocks']);
    for ii=-nlags:nlags,
        sigmaXCF(ii+nlags+1,1) = sqrt(1/(n-abs(ii)));
        BoundsXCF(ii+nlags+1,[1 2]) = sigmaXCF(ii+nlags+1,1)*[2 -2];
    end
    for j=1:length(variables_)
        % autocorr/crosscorr: 0: prefilter: assume data are without the mean,
        % consistent with moments computed around the steady state of  the
        % model
        [ACF,Lags] = autocorr(data(:,j),nlags,0,2,0);
        for ii=1:nlags,
            sigmaQ(ii,1) = sqrt((1 + 2*(ACF(2:ii)'*ACF(2:ii)))/n);
            Bounds(ii,[1 2]) = sigmaQ(ii,1)*[2 -2];
        end
        out.autocorr(j,:)=ACF(2:end);
        out.acf.(variables_{j})=[Lags(2:end) ACF(2:end) Bounds];
        %
        %  Plot the sample ACF.
        %
        %    figure(hACF)
        set(0,'CurrentFigure',hACF)
        subplot(nrow,ncol,j)
        lineHandles  =  stem([1:nlags] , ACF(2:end),'filled' , 'b-o');
        set   (lineHandles(1) , 'MarkerSize' , 2)
        xlabel('Lag')
        ylabel('ACF')
        title (variables_{j},'interpreter','none')
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
        
        plot([1:nlags] , Bounds , '--r');
        plot([1 nlags] , [0 0] , '-b');
        hold('off')
        a  =  axis;
        axis([0.5 nlags a(3) 1]);
        hCCF = dyn_figure(options_.nodisplay,'Name',[variables_{j} ' - CCF''s - shocks']);
        for i=1:length(variables_)
            [XCF,Lags] = crosscorr(data(:,j),data(:,i),nlags,2,0);
            
            %
            %  Plot the sample XCF.
            %
            subplot(nrow,ncol,i)
            lineHandles  =  stem(Lags , XCF , 'filled' , 'b-o');
            set   (lineHandles(1) , 'MarkerSize' , 2)
            xlabel('lag','fontsize',7)
            ylabel('XCF','fontsize',7)
            title ([variables_{j} '(t)  ' variables_{i} '(t-lag)'],'interpreter','none','fontsize',5)
            hold  ('on')
            %
            %  Plot the confidence bounds under the hypothesis
            %  that the underlying series are uncorrelated.
            %
            plot([-nlags:nlags] , BoundsXCF , '--r');
            
            plot([-nlags nlags] , [0 0] , '-b');
            hold('off')
            a  =  axis;
            axis([-nlags-0.5 nlags+0.5 a(3:end)]);
            set(gca,'fontsize',7)
            out.ccf.([variables_{j} '_' variables_{i}]) = [Lags XCF BoundsXCF];
        end
        set(hCCF, ...
            'position',posit, ...
            'paperorientation','landscape', ...
            'papertype','a3', ...
            'PaperPositionMode','auto')
        dyn_saveas(hCCF,[DirectoryName,filesep,fname_,'_CCF_shocks_',variables_{j}],options_.nodisplay,options_.graph_format)
    end
    set(hACF, ...
        'position',posit, ...
        'paperorientation','landscape', ...
        'papertype','a3', ...
        'PaperPositionMode','auto')
    dyn_saveas(hACF,[DirectoryName,filesep,fname_,'_ACF_shocks'],options_.nodisplay,options_.graph_format)
    
    
end

save([M_.fname '_moments_shocks.mat'],'out','names','variables_')
