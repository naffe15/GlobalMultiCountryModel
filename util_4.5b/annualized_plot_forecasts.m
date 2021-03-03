function  AnnForecastPlot = annualized_plot_forecasts(pplotvar, q2a, M_, options_, oo_, T, Tlim, nplots, pplottitles, occbin, cfrcst, frcst_ext, nfrcst_ext, ecfin_yes)


% Setting the time
t0=min(find(T==floor(T)));
forecast_ = ceil(options_.forecast/4);
T = T(t0:4:end);
TF  = T(end - forecast_ +1 : end);

if frcst_ext
    TFext = T(end)+1:T(end)+(nfrcst_ext/4);
    TFext = TFext';
    TTFext = [T; TFext];
    Tplot = TTFext;
    TF = TFext;
else
    Tplot = T;
end



% this option allows to define the forecast plot time spans
if nargin<6 || isempty(Tlim)
    % Tlim=[TT(1) TT(end)];
    Tlim=[2005 2020];
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
% this option allows to plot the forecasts with the OCCBIN toolbox.
if nargin<9 || isempty(occbin),
    occbin=0;
end
% this option allows to plot the ECFIN forecasts.
if nargin<10 || isempty(cfrcst),
    cfrcst=0;
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
    nbofplots = [nbofplots repmat(nbofplots(end),1,nfigplus)];
end

SmoothedVariables=oo_.SmoothedVariables;


bvar_yes     = 0;
occbin_yes   = 0;
exopath_yes  = 0;
endopath_yes = 0;
exog_assmpt_yes=0;
% setting the names of the figure to save

if frcst_ext ==0
    fnam_suffix = '_annual_frcst';
else
    fnam_suffix = '_annual_frcst_ext';
end

if isempty(strmatch('frcst_exog_assmpt',fieldnames(oo_.jrc)))== 0 
    fnam_suffix = [ fnam_suffix '_exog_assmpt'];
    exog_assmpt_yes=1;
end

if isempty(strmatch('bvar',fieldnames(oo_)))== 0
    fnam_suffix = [ fnam_suffix '_bvar'];
    bvar_yes=1;
end
if isempty(strmatch('forecast_exo_path',fieldnames(oo_)))==0
    fnam_suffix = [ fnam_suffix '_exopath'];
    exopath_yes=1;
end
if isempty(strmatch('forecast_endo_path',fieldnames(oo_)))==0
    fnam_suffix = [ fnam_suffix '_endopath'];
    endopath_yes=1;
end
if isempty(strmatch('occbin_forecast',fieldnames(oo_)))==0
    fnam_suffix = [ fnam_suffix '_occbin'];
    occbin_yes =1;
end
if nargin<14 || isempty(ecfin_yes),
    if cfrcst==1
        ecfin_yes = 1;
    else
        ecfin_yes = 0;
    end
end
if cfrcst==1 && ecfin_yes,
    fnam_suffix = [ fnam_suffix '_ecfin'];
end
fnam_suffix = [ fnam_suffix '_'];
% delete(['*',fnam_suffix,'*.*']);


npreamble = (forecast_+1)*4-options_.forecast;
for j=1:length(pplotvar)
    ypreamble = SmoothedVariables.(pplotvar{j})(end-(forecast_+1)*4+1:end);
    
    ypreamble = ypreamble(1:npreamble);
    % unconditional forecasts
    ytmp = [ypreamble; oo_.forecast.Mean.(pplotvar{j})]-get_mean(pplotvar{j});
    if isstruct(q2a(j).aux),
         ytmp_aux_name = q2a(j).aux.y;
         ytmp_aux_smooth(:,j) = SmoothedVariables.(ytmp_aux_name)(t0:end)-get_mean(ytmp_aux_name);
         ytmp_aux = SmoothedVariables.(ytmp_aux_name)(end-(forecast_+1)*4+1:end);
         ytmp_aux = ytmp_aux(1:npreamble);
         ytmp_aux = [ytmp_aux; oo_.forecast.Mean.(ytmp_aux_name)]-get_mean(ytmp_aux_name);
         q2a(j).aux.y = ytmp_aux;
    end
    [ya, yass, gya, gyass] = ...
        quarterly2annual(ytmp,get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
    if q2a(j).plot ==1,
        ztmp=gya(2:end)+gyass;
    elseif q2a(j).plot == 2
        ztmp=ya(2:end)+yass;
    else
        error('choose level or growth rate')
    end
    forecastedvar(:,j,1) = ztmp;
    
    if ~isstruct(q2a(j).aux),
        
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
    else
        ztmp = ztmp*nan;
    end
    forecastedvar(:,j,2) = ztmp;
    
    if ~isstruct(q2a(j).aux),        
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
    else
        ztmp = ztmp*nan;
    end
    forecastedvar(:,j,3) = ztmp;
    
    if isstruct(q2a(j).aux),
        q2a(j).aux.y = ytmp_aux_smooth(:,j);
    end
    [ya, yass, gya, gyass] = ...
        quarterly2annual(SmoothedVariables.(pplotvar{j})(t0:end)-get_mean(pplotvar{j}),get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
    if q2a(j).plot ==1,
        ztmp=gya+gyass;
    elseif q2a(j).plot == 2
        ztmp=ya+yass;
    else
        error('choose level or growth rate')
    end
    smoothedvar(:,j)=ztmp;
    
    %Model forecast with exog assumptions
    
     if exog_assmpt_yes == 1
  
       ypreamble = SmoothedVariables.(pplotvar{j})(end-(forecast_+1)*4+1:end);
    
       ypreamble = ypreamble(1:npreamble);
       ytmp = [ypreamble; oo_.jrc.frcst_exog_assmpt.SmoothedVariables.(pplotvar{j})(end-(options_.forecast-1):end)]-get_mean(pplotvar{j});
       if isstruct(q2a(j).aux),
           ytmp_aux = SmoothedVariables.(ytmp_aux_name)(end-(forecast_+1)*4+1:end);
           ytmp_aux = ytmp_aux(1:npreamble);
           ytmp_aux = [ytmp_aux; oo_.jrc.frcst_exog_assmpt.SmoothedVariables.(ytmp_aux_name)(end-(options_.forecast-1):end)]-get_mean(ytmp_aux_name);
           q2a(j).aux.y = ytmp_aux;
       end
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        
        smoothvar_exog_assmpt(:,j)=ztmp;
       
        
    end
    
    % FF added the following line
    if bvar_yes == 1 && isempty(strmatch(pplotvar{j},fieldnames(oo_.bvar.forecast.no_shock.Mean)))==0
        if isstruct(q2a(j).aux),
            q2a(j).aux.y = [ytmp_aux(1:npreamble); oo_.bvar.forecast.no_shock.Mean.(ytmp_aux_name)-get_mean(ytmp_aux_name)];
        end

        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.bvar.forecast.no_shock.Mean.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        forecastedbvar(:,j,1)=ztmp;
        
        if ~isstruct(q2a(j).aux),
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.bvar.forecast.no_shock.HPDsup.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        else
            ztmp = ztmp*nan;
        end
        forecastedbvar(:,j,2)=ztmp;
        if ~isstruct(q2a(j).aux),
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.bvar.forecast.no_shock.HPDinf.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        end
        forecastedbvar(:,j,3)=ztmp;
    end
    if exopath_yes == 1
        if isstruct(q2a(j).aux),
            q2a(j).aux.y = [ytmp_aux(1:npreamble); oo_.forecast_exo_path.Mean.(ytmp_aux_name)-get_mean(ytmp_aux_name)];
        end
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.forecast_exo_path.Mean.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        forecastedexop(:,j,1)=ztmp;
        if ~isstruct(q2a(j).aux),
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.forecast_exo_path.HPDsup.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        else
            ztmp = ztmp*nan;
        end
        forecastedexop(:,j,2)=ztmp;
        if ~isstruct(q2a(j).aux),
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.forecast_exo_path.HPDinf.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        end
        forecastedexop(:,j,3)=ztmp;
    end
    if endopath_yes == 1
        if isstruct(q2a(j).aux),
            q2a(j).aux.y = [ytmp_aux(1:npreamble); oo_.forecast_endo_path.Mean.(ytmp_aux_name)-get_mean(ytmp_aux_name)];
        end
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.forecast_endo_path.Mean.(pplotvar{j})(2:end)]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        forecastedendop(:,j,1)=ztmp;
        if ~isstruct(q2a(j).aux),
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.forecast_endo_path.ci.(pplotvar{j})(2,2:end)']-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        else
            ztmp = ztmp*nan;
        end
        forecastedendop(:,j,2)=ztmp; % sup

        if ~isstruct(q2a(j).aux),
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.forecast_endo_path.ci.(pplotvar{j})(1,2:end)']-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        end
        forecastedendop(:,j,3)=ztmp; % inf
    end
    if occbin_yes == 1
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.occbin_forecast.Mean.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        forecastedoccb(:,j,1)=ztmp;
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.occbin_forecast.HPDsup.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        forecastedoccb(:,j,2)=ztmp;
        [ya, yass, gya, gyass] = ...
            quarterly2annual([ypreamble; oo_.occbin_forecast.HPDinf.(pplotvar{j})]-get_mean(pplotvar{j}), ...
            get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
        if q2a(j).plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a(j).plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
        forecastedoccb(:,j,3)=ztmp;
    end
    if q2a(j).plot ==1,
        steadyvar(j)=gyass;
    elseif q2a(j).plot == 2
        steadyvar(j)=yass;
    else
        error('choose level or growth rate')
    end
end



jfig   = 0;
jplot  = 0;
is     = find(T==Tlim(1));

for j=1:length(pplotvar)
    xlsinfo={};
    
    varnamej=pplotvar{j};
    varpos= strmatch(varnamej, M_.endo_names, 'exact');
    texname{j,:}=q2a(j).tex_gname;
    varname{j,:}=varnamej;
    
    
    if jplot==0,
        dyn_figure(options_.nodisplay,'name','Projection plots');
        jfig=jfig+1;
    end
    jplot=jplot+1;
    subplot(nplots(jfig,1),nplots(jfig,2),jplot)
    
    if isstruct(q2a(j).aux),
        q2a(j).aux.y = ytmp_aux_smooth(:,j);
    end
    
    [ya, yass, gya, gyass] = ...
        quarterly2annual(SmoothedVariables.(pplotvar{j})(t0:end)-get_mean(pplotvar{j}),get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
    if q2a(j).plot ==1,
        ztmp=gya+gyass;
    elseif q2a(j).plot == 2
        ztmp=ya+yass;
    else
        error('choose level or growth rate')
    end
    EC_plot            = ztmp;
    actual_values      = EC_plot(1:length(T) - forecast_);
    
    if frcst_ext
        ECfrcst      = EC_plot (length(T) - forecast_+1:end);
    else
        ECfrcst = [];
    end
    
    % DSGE uncoditional forecast, mean value
    vplot           =  [actual_values ; forecastedvar(:,j,1)];
    if frcst_ext
        vplot1_           = [actual_values; NaN*(ECfrcst); NaN*(forecastedvar(:,j,1))];
        vplot2_           = [NaN*(actual_values(1:end)); ECfrcst; NaN*(forecastedvar(:,j,1))];
        vplot3_           = [NaN*(actual_values); NaN*(ECfrcst(1:end)); forecastedvar(:,j,1)];
        
        % Effectively plotted
        
        vplot1_p           = vplot1_(is:end);
        vplot2_p           = vplot2_(is:end);
        vplot3_p           = vplot3_(is:end);
        
        %find last not NaN
        ind_last1 = length(vplot1_p(~isnan(vplot1_p)));
        vplot2_p(ind_last1) = vplot1_p((ind_last1));
        ind_last2 = sum(isnan(vplot3_p));
        vplot3_p(ind_last2)= vplot2_p( ind_last2 );
        
        vplot1           = vplot1_;
        vplot2           = vplot2_;
        vplot2(length(actual_values)) = actual_values(end);
        vplot3           = vplot3_;
        vplot3(length(actual_values)+length(ECfrcst)) = ECfrcst(end);
        %vplot3           = [NaN*(actual_values); NaN*(ECfrcst(1:end-1)); ECfrcst(end); forecastedvar(:,j,1)];
    end
    
    if occbin_yes == 1
        % DSGE occbin forecast, max
        vplotStd1       = [actual_values ; ECfrcst; forecastedoccb(:,j,2)];
        % DSGE occbin forecast, min
        vplotStd_1      = [actual_values ; ECfrcst; forecastedoccb(:,j,3)];
    else
        % DSGE uncoditional forecast, upper bound
        vplotStd1       = [actual_values ; ECfrcst; forecastedvar(:,j,2)];
        % DSGE uncoditional forecast, lower bound
        vplotStd_1      = [actual_values ; ECfrcst; forecastedvar(:,j,3)];
    end
    
    % plot confidence sets
    %         floor1 = floor((length(Tplot))/4);
    %         % Get the whole years over which the quarterly data is spanning
    %         if(mod(length(Tplot), 4)>0)
    %             floor1 = floor1 +1;
    %         end
    %         %floor2 = floor((length(vplotStd_1)-3)/4);
    %         %if floor1 ~= floor2
    %             endQA = floor1*4 - mod(length(Tplot), 4);
    %         %else
    %         %    endQA = length(Tplot);
    %         %end
    TimeLineQA = [is:length(vplotStd_1)];
    %h = area([Tplot(is:4:end)],[vplotStd_1(is+3:4:end), vplotStd1(is+3:4:end) - vplotStd_1(is+3:4:end)]);
    h = area([Tplot(TimeLineQA)],[vplotStd_1(is:end), vplotStd1(is:end) - vplotStd_1(is:end)]);
    %h = area([Tplot(is:4:end-3)],[vplotStd_1(is+3:4:end), vplotStd1(is+3:4:end) - vplotStd_1(is+3:4:end)]);
    set(h(2),'FaceColor',[.95 .95 .95])
    set(h(1),'FaceColor',[1 1 1])
    hold on;
    % GM occbin forecasts
    if isfield(oo_,'MeanForecast')
        my_forecast_ = length(oo_.MeanForecast.HPDsup.(pplotvar{j}));
        ypreamble = SmoothedVariables.(pplotvar{j})(end-(ceil(my_forecast_/4)+1)*4+1:end);
        npreamble0 = (ceil(my_forecast_/4)+1)*4-my_forecast_;
    
    ypreamble0 = ypreamble(1:npreamble0);

    [ya, yass, gya, gyass] = ...
        quarterly2annual([ypreamble0; oo_.MeanForecast.HPDinf.(pplotvar{j})]-get_mean(pplotvar{j}),get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
    if q2a(j).plot ==1,
        zinf=gya(2:end)+gyass;
    elseif q2a(j).plot == 2
        zinf=ya(2:end)+yass;
    else
        error('choose level or growth rate')
    end
    [ya, yass, gya, gyass] = ...
        quarterly2annual([ypreamble0; oo_.MeanForecast.HPDsup.(pplotvar{j})]-get_mean(pplotvar{j}),get_mean(pplotvar{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
    if q2a(j).plot ==1,
        zsup=gya(2:end)+gyass;
    elseif q2a(j).plot == 2
        zsup=ya(2:end)+yass;
    else
        error('choose level or growth rate')
    end
        my_forecast_ = ceil(length(oo_.MeanForecast.HPDinf.(pplotvar{j}))/4);
        fill([Tplot(TimeLineQA(end-my_forecast_+1:end)); Tplot(TimeLineQA(end:-1:end-my_forecast_+1))],[zinf;zsup(end:-1:1)],[0.75 0.75 0.95]);
        
    end
    if exog_assmpt_yes ==1
        % occbin forecast, mean value
        exogassmptplot        = [actual_values; ECfrcst; smoothvar_exog_assmpt(:,j,1)];
        hold on, plot(Tplot(TimeLineQA),exogassmptplot(is:end),'g','Linewidth',2)
    end
    if occbin_yes ==1
        % occbin forecast, mean value
        occbinplot        = [actual_values; ECfrcst; forecastedoccb(:,j,1)];
        hold on, plot(Tplot(TimeLineQA),occbinplot(is:end),'m','Linewidth',2)
    end
    % BVAR uncoditional forecast, mean value
    if bvar_yes ==1 && isempty(strmatch(pplotvar{j},fieldnames(oo_.bvar.forecast.no_shock.Mean)))==0
        bvarplot        = [actual_values; ECfrcst; forecastedbvar(:,j,1)];
        %hold on, plot(T(is:end),bvarplot(is:end),'r','Linewidth',2)
        hold on, plot(Tplot(TimeLineQA),bvarplot(is:end),'r','Linewidth',2)
    end
    % GM forecast conditional on exogenous variables paths
    if exopath_yes ==1
        % Exogenous Paths uncoditional forecast, mean value
        exopplot        = [actual_values; ECfrcst; forecastedexop(:,j,1)];
        % hold on, plot(T(is:end),exopplot(is:end),'g','Linewidth',2)
        hold on, plot(Tplot(TimeLineQA),exopplot(is:end),'r','Linewidth',2)
    end
    % GM forecast conditional on endogenous variables paths
    if endopath_yes ==1
        % Endogenous Paths forecast, mean value
        endopplot        = [actual_values; ECfrcst; forecastedendop(:,j,1)];
        % hold on, plot(T(is:end),endopplot(is:end),'Linewidth',2,'Color',[0 0.5 0])
        hold on, plot(Tplot(TimeLineQA),endopplot(is:end),'Linewidth',2,'Color',[0 0.5 0])
    end
    % ECFIN forecast
    if frcst_ext == 0
        if ecfin_yes == 1,
            % hold on, plot(T(is:end),EC_plot(is:end),'k-.','Linewidth',2)
            hold on, plot(Tplot(TimeLineQA),EC_plot(is:end),'k-.','Linewidth',2)
        end
    end
    
    % GM unconditional forecasts
    % hold on, plot(T(is:end),vplot(is:end),'b','Linewidth',2);
    if frcst_ext == 1
        %hold on, plot(Tplot(is:4:endQA),vplot1(is+3:4:end),'b','Linewidth',2), plot(Tplot(is:4:endQA),vplot2(is+3:4:end),'k','Linewidth',2), hold on, plot(Tplot(is:4:endQA), vplot3(is+3:4:end),'k--','Linewidth',2)
        hold on, plot(Tplot(TimeLineQA),vplot1_p,'b','Linewidth',2), plot(Tplot(TimeLineQA),vplot2_p,'k','Linewidth',2), hold on, plot(Tplot(TimeLineQA), vplot3_p,'k--','Linewidth',2)
    else
        hold on, plot(Tplot(TimeLineQA),vplot(is:end),'b','Linewidth',2)
    end
    
    if q2a(j).plot ==1,
        title(q2a(j).gname,'interpreter','none')
    elseif q2a(j).plot == 2
        title(q2a(j).name,'interpreter','none')
    else
        error('choose level or growth rate')
    end
    hold on, plot([Tplot(is) Tplot(end)],[steadyvar(j) steadyvar(j)],'k:')
    axis tight
    
    %% Save the plotted figures
    if q2a(j).plot ==1
        AnnForecastPlot(j).VarName = q2a(j).gname;
    elseif q2a(j).plot == 2
        AnnForecastPlot(j).VarName = q2a(j).name;
    end    
    q2a(j).plot
    AnnForecastPlot(j).VarName
    %AnnForecastPlot(j).steadyvar    = steadyvar(j);
    AnnForecastPlot(j).TimeLineQA = Tplot(TimeLineQA);
    %AnnForecastPlot(j).TimeLine   = [Tplot(is) Tplot(end)];
    if exog_assmpt_yes ==1 % occbin forecast, mean value (green line)
        AnnForecastPlot(j).exogassmpt = exogassmptplot(is:end);
    end
    if occbin_yes ==1 % occbin forecast, mean value (magenta line)
        AnnForecastPlot(j).exogassmpt = occbinplot(is:end);
    end
    if bvar_yes ==1 && isempty(strmatch(pplotvar{j},fieldnames(oo_.bvar.forecast.no_shock.Mean)))==0    % BVAR uncoditional forecast, mean value (red line)
        AnnForecastPlot(j).bvarplot = bvarplot(is:end);
    end
    if exopath_yes ==1 % GM forecast conditional on exogenous variables paths; Exogenous Paths uncoditional forecast, mean value (red line)
        AnnForecastPlot(j).exopplot = exopplot(is:end);
    end    
    if endopath_yes ==1 % GM forecast conditional on endogenous variables paths; Endogenous Paths forecast, mean value (greenish large line)
        AnnForecastPlot(j).endopplot = endopplot(is:end);
    end
    if frcst_ext == 0 % ECFIN forecast, not extended (large black dash-dot line)
        if ecfin_yes == 1
            AnnForecastPlot(j).EC_plot = EC_plot(is:end);
        end
    end
    % GM unconditional forecasts
    if frcst_ext == 1
        AnnForecastPlot(j).Mean  = vplot1_p;
        AnnForecastPlot(j).HPD   = vplot2_p;
        AnnForecastPlot(j).EcFin = vplot3_p;
    else
        AnnForecastPlot(j).Mean = vplot(is:end);
        AnnForecastPlot(j).HPDsup = vplotStd1(is:end);
        AnnForecastPlot(j).HPDinf = vplotStd_1(is:end);
    end

        
    %%
    
    vplot0         = forecastedvar(:,j,1);
    ycomp          = forecastedvar(:,j,1)*NaN;
    yplot          = vplot0;
    yplot_sterr    = [forecastedvar(:,j,2) forecastedvar(:,j,3)];
    ysteady        = repmat(steadyvar(j),[length(TF) 1]);
    xlsinfo(1,1:6) = {'time', pplotvar{j}, '90% conf.', '90% conf.', 'steady state' , ['Comp ',pplotvar{j}]};
    xlsinfo(2:length(TF)+1,1:6) = num2cell([TF yplot yplot_sterr ysteady ycomp]);
    % xlswrite([M_.fname,'_frcst.xls'],xlsinfo,pplotvar{j});
    
    if jplot==nbofplots(jfig) || j==length(pplotvar)
        files(jfig).name = [M_.fname,fnam_suffix, int2str(Tlim(2)),'_',int2str(jfig)];
        dyn_saveas(gcf,files(jfig).name,options_.nodisplay,options_.graph_format);
        jplot=0;
    end
end


fidTeX  = fopen([M_.fname '_Forecasts.TeX'],'w');
% files   = dir(['*',fnam_suffix,'*.eps']);
h1      = 1;

numberofvars=length(pplotvar);
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

end

function [forecastedvar, ytmp_aux_smooth] = make_avar(SmoothedVariables,forecast,forecast_,npreamble,pplotvar,q2a)

    ypreamble = SmoothedVariables.(pplotvar)(end-(forecast_+1)*4+1:end);    
    ypreamble = ypreamble(1:npreamble);
    ytmp = [ypreamble; forecast.Mean.(pplotvar)]-get_mean(pplotvar);
    if isstruct(q2a.aux),
         ytmp_aux_name = q2a.aux.y;
         ytmp_aux_smooth = SmoothedVariables.(ytmp_aux_name)(t0:end)-get_mean(ytmp_aux_name);
         ytmp_aux = SmoothedVariables.(ytmp_aux_name)(end-(forecast_+1)*4+1:end);
         ytmp_aux = ytmp_aux(1:npreamble);
         ytmp_aux = [ytmp_aux; forecast.Mean.(ytmp_aux_name)]-get_mean(ytmp_aux_name);
         q2a.aux.y = ytmp_aux;
    end
    [ya, yass, gya, gyass] = ...
        quarterly2annual(ytmp,get_mean(pplotvar),q2a.GYTREND0,q2a.type,q2a.islog,q2a.aux);
    if q2a.plot ==1,
        ztmp=gya(2:end)+gyass;
    elseif q2a.plot == 2
        ztmp=ya(2:end)+yass;
    else
        error('choose level or growth rate')
    end
    forecastedvar(:,1) = ztmp;
    
    if ~isstruct(q2a.aux),
        
        ytmp = [ypreamble; forecast.HPDsup.(pplotvar)]-get_mean(pplotvar);
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar),q2a.GYTREND0,q2a.type,q2a.islog,q2a.aux);
        if q2a.plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a.plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
    else
        ztmp = ztmp*nan;
    end
    forecastedvar(:,2) = ztmp;
    
    if ~isstruct(q2a.aux),        
        ytmp = [ypreamble; forecast.HPDinf.(pplotvar)]-get_mean(pplotvar);
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar),q2a.GYTREND0,q2a.type,q2a.islog,q2a.aux);
        if q2a.plot ==1,
            ztmp=gya(2:end)+gyass;
        elseif q2a.plot == 2
            ztmp=ya(2:end)+yass;
        else
            error('choose level or growth rate')
        end
    else
        ztmp = ztmp*nan;
    end
    forecastedvar(:,3) = ztmp;

end
