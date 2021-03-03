function  [r2a, vname] = annualized_plot_fit(pplotvar, q2a, M_, options_, oo_, T, Tlim, nplots, pplottitles, pposterior)


% Setting the time
TQ = T;
% if (T(end)-floor(T(end)))~=0.75
%    T=T(T<floor(T(end)));
% end
%% WARNING: to recover r2a results prior to 2019, Jan 7th, set manually options_.presample=0 before calling this utility
% t0=min(find(T==floor(T(options_.presample+1))));
t0=min(find(T==floor(T)));
%% end WARNING
forecast = size(oo_.FilteredVariablesKStepAhead,1);
forecast_ = ceil(forecast/4);
T = T(t0:4:end);
TF  = T(end - forecast_ +1 : end);
Tplot = T;



% this option allows to define the forecast plot time spans
if nargin<7 || isempty(Tlim)
    Tlim=[T(1) T(end)];
%     Tlim=[2010 2020];
end
% this option allows to set the number of subplots
if nargin<8 || isempty(nplots)
    nplots=[3 2];
    nfigs = ceil(length(pplotvar)/6);
    nplots = repmat(nplots,nfigs,1);
end

if nargin<9 || isempty(pplottitles),
    pplottitles=pplotvar;
end

if nargin<10 || isempty(pposterior),
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

fnam_suffix = '_annual_fit';

% check pplotvar list
agname = [{q2a.gname}]';
aname = [{q2a.name}]';
qname = [{q2a.qname}]';
pplotvar0 = pplotvar;
q2a0=q2a;
q2a=q2a0(1:length(pplotvar));
jaux = 0;
for j=1:length(pplotvar)
    jvar = strmatch(pplotvar{j},M_.endo_names,'exact');
    if isempty(jvar)
        % try with the annualized name
        jj = strmatch(pplotvar{j},agname,'exact');
        if isempty(jj)
            jj = strmatch(pplotvar{j},aname,'exact');
        end   
        pplotvar{j} = q2a0(jj).qname;
    else
        if ~isequal(pplotvar{j},qname{j})
            jj = strmatch(pplotvar{j},qname,'exact');
        else
            jj = j;
        end
    end
    if isempty(jj)
        error(['Input variable name ' pplotvar{j} ' is not among annual or quarterly names in q2a list'])
    end
    q2a(j)=q2a0(jj(1));
    if length(jj) > 1
        warning(['Input variable name ' pplotvar{j} ' is ambiguous (more entries in q2a are possible)'])
    end
end
% end check

npreamble = (forecast_+1)*4-forecast;
for j=1:length(pplotvar)
    jvar = strmatch(pplotvar{j},M_.endo_names,'exact');
    if ~isstruct(q2a(j).aux)
    for t=2:length(T)
        ypreamble = oo_.UpdatedVariables.(pplotvar{j})(t0+(t-2)*4:t0+(t-1)*4-1);
        ypreamble = ypreamble(1:npreamble);
        qfrcst=diag(squeeze(oo_.FilteredVariablesKStepAhead(:,jvar,t0+(t-1)*4:t0+(t-1)*4+forecast-1)));
        ytmp = [ypreamble; qfrcst]-get_mean(pplotvar{j});
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        forecastedvar(:,t-1,j)=ztmp;

        if pposterior
        % unconditional forecasts
        ytmp = [ypreamble; oo_.forecast.Mean.(pplotvar{j})]-get_mean(pplotvar{j});
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        posterior_forecastedvar(:,j,1) = ztmp;
        
        ytmp = [ypreamble; oo_.forecast.HPDsup.(pplotvar{j})]-get_mean(pplotvar{j});
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        posterior_forecastedvar(:,j,2) = ztmp;
        
        ytmp = [ypreamble; oo_.forecast.HPDinf.(pplotvar{j})]-get_mean(pplotvar{j});
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        posterior_forecastedvar(:,j,3) = ztmp;
        end 
    end
    
    [ya, yass, gya, gyass] = ...
        quarterly2annual(oo_.UpdatedVariables.(pplotvar{j})(t0:end)-get_mean(pplotvar{j}),get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
    if q2a(j).plot ==1,
        ztmp=gya+gyass;
    elseif q2a(j).plot == 2
        ztmp=ya+yass;
    else
        error('choose level or growth rate')
    end
    updatedvar(:,j)=ztmp;
    
    % FF added the following line
    if q2a(j).plot ==1,
        steadyvar(j)=gyass;
    elseif q2a(j).plot == 2
        steadyvar(j)=yass;
    else
        error('choose level or growth rate')
    end
   for jfcast = 1:forecast_,
%        r2a(j,jfcast) = 1-sum((updatedvar(jfcast+1:end,j)-squeeze(forecastedvar(jfcast,1:end-(jfcast-1),j))').^2)/sum((updatedvar((jfcast+1):end,j)-steadyvar(j)).^2) ;
       r2a(j,jfcast) = 1-sum((updatedvar(jfcast+1:end-forecast_,j)-squeeze(forecastedvar(jfcast,1:end-forecast_-(jfcast-1),j))').^2)/sum((updatedvar((jfcast+1):end-forecast_,j)-steadyvar(j)).^2) ; 
   end
    end
end

jfig   = 0;
jplot  = 0;
is     = find(T==Tlim(1));

for j=1:length(pplotvar)
    xlsinfo={};
    
    if ~isstruct(q2a(j).aux)
    varnamej=pplotvar{j};
    varpos= strmatch(varnamej, M_.endo_names, 'exact');
    texname{j,:}=q2a(j).tex_gname;
    varname{j,:}=varnamej;
    
    
    if jplot==0,
        dyn_figure(options_.nodisplay,'name','Annual fit');
        jfig=jfig+1;
    end
    jplot=jplot+1;
    subplot(nplots(jfig,1),nplots(jfig,2),jplot)
    
    plot(T,updatedvar(:,j),'linewidth',2)
    hold on,
    for jt=1:size(forecastedvar,2)
        plot(T(jt):T(jt)+forecast_,[updatedvar(jt,j); forecastedvar(:,jt,j)],'r','linewidth',2)
    end
    if q2a(j).plot ==1,
        title(q2a(j).gname,'interpreter','none')
        vname{j,1} = q2a(j).gname;
    elseif q2a(j).plot == 2
        title(q2a(j).name,'interpreter','none')
        vname{j,1} = q2a(j).name;
    else
        error('choose level or growth rate')
    end
    
    plot(T(2:end-forecast_),squeeze(forecastedvar(1,1:end-forecast_,j))','g-.')
    if forecast_==2
    plot(T(3:end-forecast_),squeeze(forecastedvar(2,1:end-forecast_-1,j))','c-.')
    end
    plot([T(1) T(end)],[steadyvar(j) steadyvar(j)],'k:')
    hold off,
    axis tight
    set(gca,'xlim',Tlim)
    if jplot==nbofplots(jfig) || j==length(pplotvar)
        files(jfig).name = [M_.fname,fnam_suffix,'_',int2str(jfig)];
        dyn_saveas(gcf,files(jfig).name,options_.nodisplay,options_.graph_format);
        jplot=0;
    end
    else
    if j==length(pplotvar) && jplot
        files(jfig).name = [M_.fname,fnam_suffix,'_',int2str(jfig)];
        dyn_saveas(gcf,files(jfig).name,options_.nodisplay,options_.graph_format);
        jplot=0;
    end
    end
end

jfig   = 0;
jplot  = 0;
for j=1:length(pplotvar)
    xlsinfo={};
    
    if ~isstruct(q2a(j).aux)
    varnamej=pplotvar{j};
    varpos= strmatch(varnamej, M_.endo_names, 'exact');
    
    
    if jplot==0,
        dyn_figure(options_.nodisplay,'name','Recursive forecasts');
        jfig=jfig+1;
    end
    jplot=jplot+1;
    subplot(nplots(jfig,1),nplots(jfig,2),jplot)
    
    plot(TQ,oo_.UpdatedVariables.(pplotvar{j}),'linewidth',2)
    hold on,
%     for jt=1:length(oo_.UpdatedVariables.(pplotvar{j}))
    for jt=1:4:length(oo_.UpdatedVariables.(pplotvar{j}))
        qforecastedvar(:,jt,j)=diag(squeeze(oo_.FilteredVariablesKStepAhead(:,varpos,jt+1:jt+forecast)));
        plot(TQ(jt):0.25:(TQ(jt)+forecast/4),[oo_.UpdatedVariables.(pplotvar{j})(jt); qforecastedvar(:,jt,j)],'r','linewidth',2)
    end
    title(pplotvar{j},'interpreter','none')
    hold on, plot([T(1) T(end)],[get_mean(pplotvar{j}) get_mean(pplotvar{j})],'k:')
    axis tight
    set(gca,'xlim',Tlim)
    if jplot==nbofplots(jfig) || j==length(pplotvar)
        files_name = [M_.fname,'_recursive_forecasts_',int2str(jfig)];
        dyn_saveas(gcf,files_name,options_.nodisplay,options_.graph_format);
        jplot=0;
    end
    else
    if j==length(pplotvar) && jplot
        files_name = [M_.fname,'_recursive_forecasts_',int2str(jfig)];
        dyn_saveas(gcf,files_name,options_.nodisplay,options_.graph_format);
        jplot=0;
    end
    end
end

fidTeX  = fopen([M_.fname '_annual_fit.TeX'],'w');
% files   = dir(['*',fnam_suffix,'*.eps']);
h1      = 1;

numberofvars=length(varname);
for jj = 1:length(files),
    
    frcstfilename=files(jj).name;
    [a,b,c]=fileparts(frcstfilename);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    for kk=1:min(numberofvars,nbofplots(jj))
        %         hh=9*(h1-1)+1;
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(varname{h1,:}),['$' deblank(texname{h1,:}) '$']);
        h1=h1+1;
    end
    numberofvars=numberofvars-nbofplots(jj);
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,['\\includegraphics[width=0.80\\textwidth] {' b '} \n']);
    fprintf(fidTeX,'\\caption{Forecast}');
    fprintf(fidTeX,'\\label{Fig:Forecast:%s}\n',int2str(jj));
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,' \n');
end

fclose(fidTeX);


