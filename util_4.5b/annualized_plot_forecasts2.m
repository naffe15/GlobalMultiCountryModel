function  AnnForecastPlot = annualized_plot_forecasts2(pplotvar, q2a, M_, options_, oo_, T, Tlim, nplots, pplottitles, occbin, cfrcst, frcst_ext, nfrcst_ext, ecfin_yes, assmpt_range2, HPDplot2, HPDplotPar)

TINPUT = T;
TINPUTlim = Tlim;
if isempty(pplotvar)
    do_the_plot=0;
    pplotvar = [{q2a.qname}]';
else
    do_the_plot=1;
end
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

% check pplotvar list
agname = [{q2a.gname}]';
aname = [{q2a.name}]';
qname = [{q2a.qname}]';
pplotvar0 = pplotvar;
qvar = qname;
q2a0=q2a;
q2a=q2a0(1:length(pplotvar));
for j=1:length(pplotvar)
    jvar = strmatch(pplotvar{j},M_.endo_names,'exact');
    jj = strmatch(pplotvar{j},agname,'exact');
    if isempty(jj)
        jj = strmatch(pplotvar{j},aname,'exact');
    end   
    if ~isempty(jvar) && ~isempty(jj)
        jx = strmatch(pplotvar{j},qname,'exact');
        if jx~=jj 
            error(['Input variable name ' pplotvar{j} ' is ambiguous: corresponds to both Q and A names. Rename one of the two.'])
        else
            warning(['Input variable name ' pplotvar{j} ' is ambiguous: Q and A names are identical. Rename one of the two.'])
            jvar=[];
        end
    end
    if isempty(jvar) && ~isempty(jj)
        % try with the annualized name
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
    qvar{jj} = pplotvar{j};
end
% end check


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

if nargin<15 || isempty(assmpt_range2),
    assmpt_range2=0;
end
if nargin<16 || isempty(HPDplot2),
    HPDplot2=1;
end
if nargin<17 || isempty(HPDplotPar),
    HPDplotPar=1;
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

SmoothedVariables=oo_.SmoothedVariables;


bvar_yes     = 0;
occbin_yes   = occbin;
exopath_yes  = 0;
endopath_yes = 0;
exog_assmpt_yes=0;
two_steps_yes =0;
% setting the names of the figure to save

if frcst_ext ==0
    fnam_suffix = '_annual_frcst';
else
    fnam_suffix = '_annual_frcst_ext';
end

if isfield(oo_.jrc,'frcst_exog_assmpt')
    fnam_suffix = [ fnam_suffix '_exog_assmpt'];
    exog_assmpt_yes=1;
end
if isfield(oo_.jrc,'forecast_second_step')
    fnam_suffix = [ fnam_suffix '_two_steps'];
    two_steps_yes=1;
end

if isfield(oo_,'bvar')
    fnam_suffix = [ fnam_suffix '_bvar'];
    bvar_yes=1;
end
if isfield(oo_.jrc,'forecast_exo_path')
    fnam_suffix = [ fnam_suffix '_exopath'];
    exopath_yes=1;
end
if isfield(oo_.jrc,'forecast_endo_path')
    fnam_suffix = [ fnam_suffix '_endopath'];
    endopath_yes=1;
end
if isfield(oo_,'occbin_forecast')
    fnam_suffix = [ fnam_suffix '_occbin'];
    occbin_yes =1;
end
if nargin<14 || isempty(ecfin_yes)
    if cfrcst==1
        ecfin_yes = 1;
    else
        ecfin_yes = 0;
    end
end
if cfrcst==1 && ecfin_yes
    fnam_suffix = [ fnam_suffix '_ecfin'];
end

assmpt_range_yes=0;

if ~isfield(options_.jrc,'assmpt_range')
    options_.jrc.assmpt_range=0;
end
if options_.jrc.assmpt_range~=0
    assmpt_range_yes=1;
end
if assmpt_range2 == 1
    assmpt_range_yes = 2;
end

% fnam_suffix = [ fnam_suffix '_'];
% delete(['*',fnam_suffix,'*.*']);


npreamble = (forecast_+1)*4-options_.forecast;
for j=1:length(pplotvar)
    [forecastedvar(:,j,:), smoothedvar(:,j), steadyvar(j), tmp] = make_avar(oo_.forecast,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),1,HPDplot2);    % unconditional forecasts
    if isstruct(q2a(j).aux),
        ytmp_aux_smooth(:,j) = tmp;
        ytmp_aux_name = q2a(j).aux.y;
    end

    if assmpt_range_yes == 1
        assmpt_range=options_.jrc.assmpt_range;
        for jjj=1:size(M_.jrc.assmpt_range_vars,1)
            
            for jj=1+(jjj-1)*assmpt_range:jjj*assmpt_range
                smoothvar_range_assmpt(:,j,jj) = make_avar(oo_.jrc.frcst_exog_assmpt_range.( M_.jrc.cfrcst_range_name{jj} ).SmoothedVariables,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),1,HPDplot2);    % unconditional forecasts
            end
            
        end
    end
    
    if assmpt_range_yes == 2
        assmpt_range=options_.jrc.assmpt_range;
        %for jjj=1:size(M_.jrc.assmpt_range_vars,1)
        jj = 0;
             for j0 = 1:length(M_.jrc.assmpt_range_vars)
                for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
                    if j2 ~= j0 
                        for j3 = 1:assmpt_range
                           for j4 = 1:assmpt_range
                               jj = jj+1;
                                 smoothvar_range_assmpt(:,j,jj)= make_avar(oo_.jrc.frcst_exog_assmpt_range.( M_.jrc.cfrcst_range_name{jj} ).SmoothedVariables,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),1,HPDplot2);
                           end
                        end
                    end
                end
             end
       
    end
    
    %Model forecast with exog assumptions
    
    if exog_assmpt_yes == 1
        smoothvar_exog_assmpt(:,j) = make_avar(oo_.jrc.frcst_exog_assmpt.SmoothedVariables,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),0);    % unconditional forecasts
    end
    if two_steps_yes == 1
        smoothvar_two_steps(:,j) = make_avar(oo_.jrc.forecast_second_step.SmoothedVariables,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),0);    % unconditional forecasts
    end
    % FF added the following line
    if bvar_yes == 1 && isempty(strmatch(pplotvar{j},fieldnames(oo_.bvar.forecast.no_shock.Mean)))==0
        forecastedbvar(:,j,:) = make_avar(oo_.bvar.forecast.no_shock,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),1,HPDplot2);    % unconditional forecasts
    end
    if exopath_yes == 1
        forecastedexop(:,j,:) = make_avar(oo_.jrc.forecast_exo_path,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),1,HPDplot2);    % unconditional forecasts
    end
    if endopath_yes == 1
	%modified by Olga to trap the nfrcst observations needed for make_avar function (and not nfrcst+1)
        frcstvar=oo_.jrc.forecast_endo_path.Mean;
        frcstvar.ci=oo_.jrc.forecast_endo_path.ci; % sup
        forecastedendop(:,j,:) = make_avar(frcstvar,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),1,HPDplot2);    % unconditional forecasts
    end
    if occbin_yes == 1
        if length(oo_.occbin.forecast.HPDinf.(pplotvar{1})) == length(oo_.occbin.forecast.Mean.(pplotvar{1}))
            HPDoccbin = 1;
        else
            HPDoccbin=0;
        end
        forecastedoccb(:,j,:) = make_avar(oo_.occbin.forecast,SmoothedVariables,options_.forecast,forecast_,npreamble,t0,pplotvar{j},q2a(j),1,HPDoccbin);    % unconditional forecasts
    end
    
end


is     = find(T==Tlim(1));

if  isfield(M_.jrc,'assmpt_range_vars')
    if assmpt_range_yes == 1
        nfigs=size(M_.jrc.assmpt_range_vars,1);
    elseif assmpt_range_yes == 2
        vvv=length(M_.jrc.assmpt_range_vars);
        nfigs=nchoosek(vvv,2);
    end
else
    nfigs=1;
end

fignam_suffix='';
for jjj=1:nfigs
    if assmpt_range_yes == 1
        fnam_suffix = [ '_annual_frcst_assmpt_range_' M_.jrc.assmpt_range_vars{jjj}];
        fignam_suffix = [': assumption range ' M_.jrc.assmpt_range_vars{jjj}];
    end
    
    if assmpt_range_yes == 2
            if jjj<4
                jj0=1;
                jj2=jj0+jjj;
                %for j0 = 1:length(M_.jrc.assmpt_range_vars)-1
                %for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
                fnam_suffix = [ '_annual_frcst_assmpt_range_' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
                fignam_suffix = [': assumption range ' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
            else
               if jjj<6
                jj0=2;
                jj2=jjj-1;
                %for j0 = 1:length(M_.jrc.assmpt_range_vars)-1
                %for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
                fnam_suffix = [ '_annual_frcst_assmpt_range_' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
                fignam_suffix = [': assumption range ' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
               else
                jj0=3;
                jj2=jj0+1;
                %for j0 = 1:length(M_.jrc.assmpt_range_vars)-1
                %for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
                fnam_suffix = [ '_annual_frcst_assmpt_range_' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
                fignam_suffix = [': assumption range ' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
               end
            end
    end
        
jfig   = 0;
jplot  = 0;
for j=1:length(pplotvar)
    xlsinfo={};
    
    varnamej=pplotvar{j};
    varpos= strmatch(varnamej, M_.endo_names, 'exact');
    texname{j,:}=q2a(j).tex_gname;
    varname{j,:}=varnamej;
    
    
    if jplot==0 && do_the_plot
        dyn_figure(options_.nodisplay,'name',['Projection plots' fignam_suffix]);
        jfig=jfig+1;
    end
    jplot=jplot+1;
    if do_the_plot
        subplot(nplots(jfig,1),nplots(jfig,2),jplot)
    end
    
    EC_plot            = smoothedvar(:,j);
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
    %h = area([Tplot(is:4:end-3)],[vplotStd_1(is+3:4:end), vplotStd1(is+3:4:end) - vplotStd_1(is+3:4:end)]);
%     h = area([Tplot(TimeLineQA)],[vplotStd_1(is:end), vplotStd1(is:end) - vplotStd_1(is:end)]);
%     set(h(2),'FaceColor',[.95 .95 .95])
%     set(h(1),'FaceColor',[1 1 1])
    if do_the_plot
    fill([Tplot(TimeLineQA); Tplot(TimeLineQA(end:-1:1))],[vplotStd_1(is:end); vplotStd1(end:-1:is)],[0.95 0.95 0.95]);
    hold on;
    
    % GM occbin forecasts
    if isfield(oo_,'MeanForecast')
        my_qforecast_ = length(oo_.MeanForecast.HPDsup.(pplotvar{j}));
        npreamble0 = (ceil(my_qforecast_/4)+1)*4-my_qforecast_;
        my_forecast_ = ceil(my_qforecast_/4);
        my_forecastedvar = make_avar(oo_.MeanForecast,SmoothedVariables,my_qforecast_,my_forecast_,npreamble0,t0,pplotvar{j},q2a(j),1,HPDplotPar);    % unconditional forecasts
        zsup = my_forecastedvar(:,2);
        zinf = my_forecastedvar(:,3);
        fill([Tplot(TimeLineQA(end-my_forecast_+1:end)); Tplot(TimeLineQA(end:-1:end-my_forecast_+1))],[zinf;zsup(end:-1:1)],[0.75 0.75 0.95]);
        
    end
    end
    if exog_assmpt_yes ==1
        % occbin forecast, mean value
        exogassmptplot        = [actual_values; ECfrcst; smoothvar_exog_assmpt(:,j,1)];
        if do_the_plot
            hold on, plot(Tplot(TimeLineQA),exogassmptplot(is:end),'g','Linewidth',2)
        end
    end
    if two_steps_yes ==1
        % occbin forecast, mean value
        twostepsplot        = [actual_values; ECfrcst; smoothvar_two_steps(:,j,1)];
        if do_the_plot
            hold on, plot(Tplot(TimeLineQA),twostepsplot(is:end),'y','Linewidth',2)
        end
    end
    if occbin_yes ==1
        % occbin forecast, mean value
        occbinplot        = [actual_values; ECfrcst; forecastedoccb(:,j,1)];
        if do_the_plot        
            hold on, plot(Tplot(TimeLineQA),occbinplot(is:end),'m','Linewidth',2)
        end
    end
    % BVAR uncoditional forecast, mean value
    if bvar_yes ==1 && isempty(strmatch(pplotvar{j},fieldnames(oo_.bvar.forecast.no_shock.Mean)))==0
        bvarplot        = [actual_values; ECfrcst; forecastedbvar(:,j,1)];
        if do_the_plot
        %hold on, plot(T(is:end),bvarplot(is:end),'r','Linewidth',2)
        hold on, plot(Tplot(TimeLineQA),bvarplot(is:end),'r','Linewidth',2)
        end
    end
    % GM forecast conditional on exogenous variables paths
    if exopath_yes ==1
        % Exogenous Paths uncoditional forecast, mean value
        exopplot        = [actual_values; ECfrcst; forecastedexop(:,j,1)];
        if do_the_plot
            % hold on, plot(T(is:end),exopplot(is:end),'g','Linewidth',2)
            hold on, plot(Tplot(TimeLineQA),exopplot(is:end),'r','Linewidth',2)
        end
    end
    % GM forecast conditional on endogenous variables paths
    if endopath_yes ==1
        % Endogenous Paths forecast, mean value
        endopplot        = [actual_values; ECfrcst; forecastedendop(:,j,1)];
        if do_the_plot
            % hold on, plot(T(is:end),endopplot(is:end),'Linewidth',2,'Color',[0 0.5 0])
            hold on, plot(Tplot(TimeLineQA),endopplot(is:end),'Linewidth',2,'Color',[0 0.5 0])
        end
    end
    % ECFIN forecast
    if frcst_ext == 0
        if ecfin_yes == 1 && assmpt_range_yes ~= 2 && do_the_plot
            % hold on, plot(T(is:end),EC_plot(is:end),'k-.','Linewidth',2)
            hold on, plot(Tplot(TimeLineQA),EC_plot(is:end),'k-.','Linewidth',2)
        end
    end
    
    if assmpt_range_yes == 1
       
       for jj=1+(jjj-1)*assmpt_range:jjj*assmpt_range
            tmp1 = [actual_values; ECfrcst];
            tmp = [tmp1(end); smoothvar_range_assmpt(:,j,jj)];
            assmpt_range_plot = [nan(length(TimeLineQA)-length(tmp),1); tmp];
            if assmpt_range == 2
                if rem(jj,2)~=0
                    if do_the_plot
                        hold on, plot(Tplot(TimeLineQA), assmpt_range_plot,'--+r','Linewidth',1)
                    end
                    AnnForecastPlot(j).assmptrangeplot.([M_.jrc.assmpt_range_vars{jjj} '_High']) = assmpt_range_plot;
                else
                    if do_the_plot
                        hold on, plot(Tplot(TimeLineQA), assmpt_range_plot,'--r','Linewidth',1)
                    end
                    AnnForecastPlot(j).assmptrangeplot.([M_.jrc.assmpt_range_vars{jjj} '_Low']) = assmpt_range_plot;
                end  
            else
                if do_the_plot
                    hold on, plot(Tplot(TimeLineQA), assmpt_range_plot,'--r','Linewidth',1)
                end
                AnnForecastPlot(j).assmptrangeplot.([M_.jrc.assmpt_range_vars{jjj} '_AdHoc']) = assmpt_range_plot;
            end
            
        end
        
    end
    
    if assmpt_range_yes == 2
        
        %  M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2}        
        
        
        for jj=1+(jjj-1)*length(M_.jrc.assmpt_range_vars):jjj*length(M_.jrc.assmpt_range_vars)
            %for j0 = 1:length(M_.jrc.assmpt_range_vars)-1
            %for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
            tmp2 = [actual_values; ECfrcst];
            tmp22 = [tmp2(end); smoothvar_range_assmpt(:,j,jj)];
            %tmp22 = smoothvar_range_assmpt(:,j,jj);
            assmpt_range_plot = [nan(length(TimeLineQA)-length(tmp22),1); tmp22];
            %assmpt_range_plot = smoothvar_range_assmpt(:,j,jj);
            if jj == 1+(jjj-1)*length(M_.jrc.assmpt_range_vars)
                %tmp_name = 'HIGH_HIGH';
                if do_the_plot
                    hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'-r','Linewidth',2)
                end
                AnnForecastPlot(j).assmptrangeplot.([M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} '_High_High']) = assmpt_range_plot;
            elseif jj == 1+(jjj-1)*length(M_.jrc.assmpt_range_vars)+1
                %tmp_name = 'HIGH_LOW';
                if do_the_plot
                    hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'--m','Linewidth',2)
                end
                AnnForecastPlot(j).assmptrangeplot.([M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} '_High_Low']) = assmpt_range_plot;
            elseif jj == 1+(jjj-1)*length(M_.jrc.assmpt_range_vars)+2
                %tmp_name = 'LOW_HIGH';
                if do_the_plot
                    hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),':c','Linewidth',2)
                end
                AnnForecastPlot(j).assmptrangeplot.([M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} '_Low_High']) = assmpt_range_plot;
            elseif jj == 1+(jjj-1)*length(M_.jrc.assmpt_range_vars)+3
                %tmp_name = 'LOW_LOW';
                if do_the_plot
                    hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'-.k','Linewidth',2)
                end
                AnnForecastPlot(j).assmptrangeplot.([M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} '_Low_Low']) = assmpt_range_plot;
            end
        end
        
    end
    
    
    % GM unconditional forecasts
    % hold on, plot(T(is:end),vplot(is:end),'b','Linewidth',2);
    if do_the_plot
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
        axis tight
        aa=axis;
        if steadyvar(j)>=aa(3)-(aa(4)-aa(3))/5 && steadyvar(j)<= aa(4)+(aa(4)-aa(3))/5
            hold on, plot([Tplot(is) Tplot(end)],[steadyvar(j) steadyvar(j)],'k:')
        end
    end
    %% Save the plotted figures
    if q2a(j).plot ==1
        AnnForecastPlot(j).VarName = q2a(j).gname; 
    elseif q2a(j).plot == 2
        AnnForecastPlot(j).VarName = q2a(j).name;  
    end
    AnnForecastPlot(j).LongName = q2a(j).frcst_name;
    %AnnForecastPlot(j).steadyvar    = steadyvar(j);
    AnnForecastPlot(j).TimeLineQA = Tplot(TimeLineQA);
    %AnnForecastPlot(j).TimeLine   = [Tplot(is) Tplot(end)];
    if exog_assmpt_yes ==1 % occbin forecast, mean value (green line)
        AnnForecastPlot(j).exogassmpt = exogassmptplot(is:end);
    end
    if two_steps_yes ==1 % occbin forecast, mean value (green line)
        AnnForecastPlot(j).twosteps = twostepsplot(is:end);
    end
    if occbin_yes ==1 % occbin forecast, mean value (magenta line)
        AnnForecastPlot(j).occbin = occbinplot(is:end);
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
    
    if do_the_plot && ((jplot==nbofplots(jfig) || j==length(pplotvar)) ) 
        if assmpt_range_yes == 2
                annotation('textbox', [0.05,0,0.2,0.05],'String', 'High High','Color','r','horizontalalignment','center','interpreter','none');
                annotation('textbox', [0.3,0,0.2,0.05],'String', 'High Low','Color','m','horizontalalignment','center','interpreter','none');
                annotation('textbox', [0.55,0,0.2,0.05],'String', 'Low High','Color','c','horizontalalignment','center','interpreter','none');
                annotation('textbox', [0.8,0,0.2,0.05],'String', 'Low Low','Color','k','horizontalalignment','center','interpreter','none');                
        end
        files(jfig).name = [M_.fname fnam_suffix '_' int2str(jfig)];
        dyn_saveas(gcf,files(jfig).name,options_.nodisplay,options_.graph_format);
        jplot=0;
    end
end
end

if do_the_plot
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

if length(pplotvar)<length(qname)
    
    AnnForecastPlot = annualized_plot_forecasts2([], q2a0, M_, options_, oo_, TINPUT, TINPUTlim, [], [], occbin, cfrcst, frcst_ext, nfrcst_ext, ecfin_yes, assmpt_range2, HPDplot2);

end

end

function [forecastedvar, smoothedvar, steadyvar, ytmp_aux_smooth] = make_avar(forecast,SmoothedVariables,qforecast_,forecast_,npreamble,t0,pplotvar,q2a,HPDplot,HPDplot2)

if nargin < 9 || isempty(HPDplot)
    HPDplot = 1;
end
if nargin < 10 || isempty(HPDplot2)
    HPDplot2 = 1;
end

ytmp_aux_smooth=[];
ypreamble = SmoothedVariables.(pplotvar)(end-(forecast_+1)*4+1:end);
ypreamble = ypreamble(1:npreamble);
if isfield(forecast,'Mean') && HPDplot
    ytmp = [ypreamble; forecast.Mean.(pplotvar)]-get_mean(pplotvar);
else
    ytmp = [ypreamble; forecast.(pplotvar)(end-(qforecast_-1):end)]-get_mean(pplotvar);
end
if isstruct(q2a.aux),
    ytmp_aux_name = q2a.aux.y;
    ytmp_aux_smooth = SmoothedVariables.(ytmp_aux_name)(t0:end)-get_mean(ytmp_aux_name);
    ytmp_aux = SmoothedVariables.(ytmp_aux_name)(end-(forecast_+1)*4+1:end);
    ytmp_aux = ytmp_aux(1:npreamble);
    if isfield(forecast,'Mean') && HPDplot
        ytmp_aux = [ytmp_aux; forecast.Mean.(ytmp_aux_name)]-get_mean(ytmp_aux_name);
    else
        ytmp_aux = [ytmp_aux; forecast.(ytmp_aux_name)(end-(qforecast_-1):end)]-get_mean(ytmp_aux_name);
    end
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
if isfield(forecast,'Mean') && HPDplot
    forecastedvar(:,2) = ztmp*nan;  
    forecastedvar(:,3) = ztmp*nan;
end
if isfield(forecast,'Mean') && HPDplot2
    if ~isstruct(q2a.aux),
        ylevel = ya;
        gylevel = gya;
        if isfield(forecast,'ci'),
            ytmp = [ypreamble; forecast.ci.(pplotvar)(2,2:end)'-get_mean(pplotvar)];
        else
            ytmp = [ypreamble; forecast.HPDsup.(pplotvar)]-get_mean(pplotvar);
        end
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar),q2a.GYTREND0,q2a.type,q2a.islog,q2a.aux);
        if q2a.plot ==1,
            % this is an uncertainty!!!!
            if q2a.islog
                ydiff = ya-ylevel;
            else
                ydiff = log(ya+yass)-log(ylevel+yass);
            end
            gya(2:end) = sqrt(ydiff(2:end).^2+ydiff(1:end-1).^2);
            ztmp=gya(2:end)+gylevel(2:end)+gyass;
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
        if isfield(forecast,'ci'),
            ytmp = [ypreamble; forecast.ci.(pplotvar)(1,2:end)'-get_mean(pplotvar)];
        else
            ytmp = [ypreamble; forecast.HPDinf.(pplotvar)]-get_mean(pplotvar);
        end
        [ya, yass, gya, gyass] = ...
            quarterly2annual(ytmp,get_mean(pplotvar),q2a.GYTREND0,q2a.type,q2a.islog,q2a.aux);
        if q2a.plot ==1,
            % this is an uncertainty!!!!
            if q2a.islog
                ydiff = ya-ylevel;
            else
                ydiff = log(ya+yass)-log(ylevel+yass);
            end
            gya(2:end) = sqrt(ydiff(2:end).^2+ydiff(1:end-1).^2);
            ztmp=-gya(2:end)+gylevel(2:end)+gyass;
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

if isstruct(q2a.aux),
    q2a.aux.y = ytmp_aux_smooth;
end
[ya, yass, gya, gyass] = ...
    quarterly2annual(SmoothedVariables.(pplotvar)(t0:end)-get_mean(pplotvar),get_mean(pplotvar),q2a.GYTREND0,q2a.type,q2a.islog,q2a.aux);
if q2a.plot ==1,
    ztmp=gya+gyass;
    steadyvar=gyass;
elseif q2a.plot == 2
    ztmp=ya+yass;
    steadyvar=yass;
else
    error('choose level or growth rate')
end
smoothedvar=ztmp;

end
