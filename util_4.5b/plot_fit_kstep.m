function  r2k = plot_fit_kstep(pplotvar, M_, options_, oo_, T, Tlim, nplots, pplottitles, pposterior)


forecast = size(oo_.FilteredVariablesKStepAhead,1);
TF  = T(end - forecast +1 : end);
Tplot = T;

r2k=[];

% this option allows to define the forecast plot time spans
if nargin<6 || isempty(Tlim)
    Tlim=[T(1) T(end)];
%     Tlim=[2010 2020];
end
% this option allows to set the number of subplots
if nargin<7 || isempty(nplots)
    nplots=[3 2];
    nfigs = ceil(length(pplotvar)/6);
    nplots = repmat(nplots,nfigs,1);
end

if nargin<8 || isempty(pplottitles),
    pplottitles=pplotvar;
end

if nargin<9 || isempty(pposterior),
    pposterior=0;
end

if size(nplots,2)==3,
    nbofplots = nplots(:,end);
    nplots=nplots(:,1:2);
else
    for j=1:size(nplots,1),
        nbofplots(j)=nplots(j,1)*nplots(j,2);
    end
end
ntotplots = sum(nbofplots);
if ntotplots<length(pplotvar),
    nfigplus = ceil((length(pplotvar)-ntotplots)/nbofplots(end));
    lastrow=nplots(end,:);
    lastrow=repmat(lastrow,nfigplus,1);
    nplots = [nplots;lastrow];
    nbofplots = [nbofplots; repmat(nbofplots(end),nfigplus,1)];
end


% setting the names of the figure to save

fnam_suffix = '_kstep_fit';

% check pplotvar list


jfig   = 0;
jplot  = 0;
for j=1:length(pplotvar)
    xlsinfo={};
    
    varnamej=pplotvar{j};
    varpos= strmatch(varnamej, M_.endo_names, 'exact');
    
    
    if jplot==0,
        dyn_figure(options_.nodisplay,'name','Recursive forecasts');
        jfig=jfig+1;
    end
    jplot=jplot+1;
    subplot(nplots(jfig,1),nplots(jfig,2),jplot)
    
    plot(T,oo_.UpdatedVariables.(pplotvar{j}),'linewidth',2)
    hold on,
%     for jt=1:length(oo_.UpdatedVariables.(pplotvar{j}))
    for jt=1:length(oo_.UpdatedVariables.(pplotvar{j}))
        if pposterior
            post_forecastedvar(1:2,1)=oo_.UpdatedVariables.(pplotvar{j})(jt);
            for jfor = 1:forecast
                post_forecastedvar(1,jfor+1)=oo_.(['Filtered_Variables_' int2str(jfor) '_step_ahead']).HPDinf.(pplotvar{j})(jt);
                post_forecastedvar(2,jfor+1)=oo_.(['Filtered_Variables_' int2str(jfor) '_step_ahead']).HPDsup.(pplotvar{j})(jt);
            end
            hf = fill([(T(jt):(T(jt)+forecast))';((T(jt)+forecast):-1:T(jt))'],[post_forecastedvar(1,:)';post_forecastedvar(2,end:-1:1)'],[0.95 0.75 0.75]); 
            set(hf,'FaceAlpha',0.3)
        end
        qforecastedvar(:,jt,j)=diag(squeeze(oo_.FilteredVariablesKStepAhead(:,varpos,jt+1:jt+forecast)));
        plot(T(jt):(T(jt)+forecast),[oo_.UpdatedVariables.(pplotvar{j})(jt); qforecastedvar(:,jt,j)],'r','linewidth',2)
    end
    for jf=1:forecast
        r2k(j,jf) = 1-sum((oo_.UpdatedVariables.(pplotvar{j})(jf+1:end)-transpose(qforecastedvar(jf,1:end-jf,j))).^2)/sum((oo_.UpdatedVariables.(pplotvar{j})(jf+1:end)-get_mean(pplotvar{j})).^2);
    end
    title(pplotvar{j},'interpreter','none')
    hold on, plot([T(1) T(end)],[get_mean(pplotvar{j}) get_mean(pplotvar{j})],'k:')
    axis tight
    set(gca,'xlim',Tlim)
    if jplot==nbofplots(jfig) || j==length(pplotvar)
        files_name = [M_.fname,fnam_suffix,int2str(jfig)];
        dyn_saveas(gcf,files_name,options_.nodisplay,options_.graph_format);
        jplot=0;
    end
end


