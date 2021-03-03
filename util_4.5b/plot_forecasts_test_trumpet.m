function  plot_forecasts_(pplotvar, M_, options_, oo_, T, Tlim, nplots, pplottitles, occbin, cfrcst, frcst_ext, nfrcst_ext, ecfin_yes, assmpt_range2)

save pplotvar pplotvar T;

% Setting the time 
TF  = T(end - options_.forecast +1 : end);

if frcst_ext
    TFext = T(end)+0.25:0.25:T(end)+(nfrcst_ext/4)';
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

if nargin<14 || isempty(assmpt_range2),
    assmpt_range2=0;
end


for j=1:size(nplots,1),
    nbofplots(j)=nplots(j,1)*nplots(j,2);
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
    fnam_suffix = '_frcst';
else
    fnam_suffix = '_frcst_ext';
end

if isfield(oo_.jrc,'frcst_exog_assmpt')
    fnam_suffix = [ fnam_suffix '_exog_assmpt'];
    exog_assmpt_yes=1;
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
if nargin<13 || isempty(ecfin_yes),
    if cfrcst==1
        ecfin_yes = 1;
    else
        ecfin_yes = 0;
    end
end
if cfrcst==1 && ecfin_yes,
    fnam_suffix = [ fnam_suffix '_ecfin'];
end

assmpt_range_yes=0;

if ~isfield(options_.jrc,'assmpt_range'),
    options_.jrc.assmpt_range=0;
end
if options_.jrc.assmpt_range~=0
    assmpt_range_yes=1;
end
if assmpt_range2 == 1
    assmpt_range_yes = 2;
end

%assmpt_range_yes = 2;
% fnam_suffix = [ fnam_suffix '_']; 
% delete(['*',fnam_suffix,'*.*']);


for j=1:length(pplotvar)
    % unconditional forecasts 
    forecastedvar(:,j,1)=oo_.forecast.Mean.(pplotvar{j});
    forecastedvar(:,j,2)=oo_.forecast.HPDsup.(pplotvar{j});
    forecastedvar(:,j,3)=oo_.forecast.HPDinf.(pplotvar{j});
    smoothedvar(:,j)=SmoothedVariables.(pplotvar{j});
    
    if assmpt_range_yes == 1
        assmpt_range=options_.jrc.assmpt_range;
        for jjj=1:size(M_.jrc.assmpt_range_vars,1)
            
            for jj=1+(jjj-1)*assmpt_range:jjj*assmpt_range
                %             eval(['load gemc_results_assmpt_',int2str(jj),' oo_ '])
                smoothvar_range_assmpt(:,j,jj)=oo_.jrc.frcst_exog_assmpt_range.( M_.jrc.cfrcst_range_name{jj} ).SmoothedVariables.(pplotvar{j});
                %             smoothvar_range_assmpt(:,j,jj)=eval(['oo_.SmoothedVariables.',pplotvar{j},';']);
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
                                 smoothvar_range_assmpt(:,j,jj)=oo_.jrc.frcst_exog_assmpt_range.( M_.jrc.cfrcst_range_name{jj} ).SmoothedVariables.(pplotvar{j});
                           end
                        end
                    end
                end
             end
       
    end
    
    
    if exog_assmpt_yes == 1
         smoothvar_exog_assmpt(:,j)=oo_.jrc.frcst_exog_assmpt.SmoothedVariables.(pplotvar{j});
    end    
    
    
    % FF added the following line
    if bvar_yes == 1 && isempty(strmatch(pplotvar{j},fieldnames(oo_.bvar.forecast.no_shock.Mean)))==0
        forecastedbvar(:,j,1)=oo_.bvar.forecast.no_shock.Mean.(pplotvar{j});
        forecastedbvar(:,j,2)=oo_.bvar.forecast.no_shock.HPDsup.(pplotvar{j});
        forecastedbvar(:,j,3)=oo_.bvar.forecast.no_shock.HPDinf.(pplotvar{j});
    end    
    if exopath_yes == 1
        forecastedexop(:,j,1)=oo_.jrc.forecast_exo_path.Mean.(pplotvar{j});
        forecastedexop(:,j,2)=oo_.jrc.forecast_exo_path.HPDsup.(pplotvar{j});
        forecastedexop(:,j,3)=oo_.jrc.forecast_exo_path.HPDinf.(pplotvar{j});
    end
    if endopath_yes == 1
        forecastedendop(:,j,1)=oo_.jrc.forecast_endo_path.Mean.(pplotvar{j})(2:end);
        forecastedendop(:,j,2)=oo_.jrc.forecast_endo_path.ci.(pplotvar{j})(2:end,2); % sup
        forecastedendop(:,j,3)=oo_.jrc.forecast_endo_path.ci.(pplotvar{j})(2:end,1); % inf
    end
    if occbin_yes == 1
        forecastedoccb(:,j,1)=oo_.occbin_forecast.Mean.(pplotvar{j});
        forecastedoccb(:,j,2)=oo_.occbin_forecast.HPDsup.(pplotvar{j});
        forecastedoccb(:,j,3)=oo_.occbin_forecast.HPDinf.(pplotvar{j});
    end
    steadyvar(j)=get_mean(pplotvar{j});
end



jfig   = 0;
jplot  = 0;
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
        if nfigs>1
            if assmpt_range_yes ==1
                fnam_suffix = [ '_frcst_assmpt_range_' M_.jrc.assmpt_range_vars{jjj}];
                fignam_suffix = [': assumption range ' M_.jrc.assmpt_range_vars{jjj}];
            end
        end
        
        if assmpt_range_yes == 2
            if jjj<4
                jj0=1;
                jj2=jj0+jjj;
                %for j0 = 1:length(M_.jrc.assmpt_range_vars)-1
                %for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
                fnam_suffix = [ '_frcst_assmpt_range_' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
                fignam_suffix = [': assumption range ' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
            else
               if jjj<6
                jj0=2;
                jj2=jjj-1;
                %for j0 = 1:length(M_.jrc.assmpt_range_vars)-1
                %for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
                fnam_suffix = [ '_frcst_assmpt_range_' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
                fignam_suffix = [': assumption range ' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
               else
                jj0=3;
                jj2=jj0+1;
                %for j0 = 1:length(M_.jrc.assmpt_range_vars)-1
                %for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
                fnam_suffix = [ '_frcst_assmpt_range_' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
                fignam_suffix = [': assumption range ' M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2} ];
               end
            end
        end
           


for j=1:length(pplotvar)
    xlsinfo={};
	Forecasts_Plots(j).Name = pplotvar{j};

    varnamej=pplotvar{j};
    varpos= strmatch(varnamej, M_.endo_names, 'exact');
    texname{j,:}=M_.endo_names_tex(varpos,:);
    varname{j,:}=varnamej;
    
    
    if jplot==0,
        dyn_figure(options_.nodisplay,'name',['Projection plots' fignam_suffix]);
        jfig=jfig+1;
    end
    jplot=jplot+1;
    subplot(nplots(jfig,1),nplots(jfig,2),jplot)
    
    EC_plot            = get_smooth(pplotvar{j});
    actual_values      = EC_plot(1:length(T) - options_.forecast);
    
    if frcst_ext
    ECfrcst      = EC_plot (length(T) - options_.forecast+1:end);    
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
    
    vplot1_p           = vplot1_(is+3:4:end);
    vplot2_p           = vplot2_(is+3:4:end);
    vplot3_p           = vplot3_(is+3:4:end);
    
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
    if findstr('QA',varnamej)~=0
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
        TimeLineQA = [is+3:4:length(vplotStd_1)]-3;
        TimeLineQ = [is+3:length(vplotStd_1)]-3;
        
        %h = area([Tplot(is:4:end)],[vplotStd_1(is+3:4:end), vplotStd1(is+3:4:end) - vplotStd_1(is+3:4:end)]);
        h = area([Tplot(TimeLineQA)],[vplotStd_1(is+3:4:end), vplotStd1(is+3:4:end) - vplotStd_1(is+3:4:end)]);
        %h = area([Tplot(is:4:end-3)],[vplotStd_1(is+3:4:end), vplotStd1(is+3:4:end) - vplotStd_1(is+3:4:end)]);
        set(h(2),'FaceColor',[.95 .95 .95])
        set(h(1),'FaceColor',[1 1 1])
        hold on;    
    else
        h = area([Tplot(is:end)],[vplotStd_1(is:end), vplotStd1(is:end) - vplotStd_1(is:end)]);
        set(h(2),'FaceColor',[.95 .95 .95])
        set(h(1),'FaceColor',[1 1 1])
        hold on;    
    end
    if isfield(oo_,'MeanForecast')
        my_forecast_ = length(oo_.MeanForecast.HPDinf.(pplotvar{j}));
        fill([Tplot(end-my_forecast_+1:end); Tplot(end:-1:end-my_forecast_+1)],[oo_.MeanForecast.HPDinf.(pplotvar{j});oo_.MeanForecast.HPDsup.(pplotvar{j})(end:-1:1)],[0.75 0.75 0.95]);
        
    end
    forecast_= options_.forecast;      
    
    % GM occbin forecasts
    if occbin_yes ==1 
        % occbin forecast, mean value
        occbinplot        = [actual_values; ECfrcst; forecastedoccb(:,j,1)];
        if findstr('QA',varnamej)~=0 
            hold on, plot(Tplot(TimeLineQA),occbinplot(is+3:4:end),'m','Linewidth',2) 
        else
            hold on, plot(Tplot(is:end),occbinplot(is:end),'m','Linewidth',2) 
        end
    end
     % BVAR uncoditional forecast, mean value
    if bvar_yes ==1 && isempty(strmatch(pplotvar{j},fieldnames(oo_.bvar.forecast.no_shock.Mean)))==0       
        bvarplot        = [actual_values; ECfrcst; forecastedbvar(:,j,1)];
        %hold on, plot(T(is:end),bvarplot(is:end),'r','Linewidth',2) 
        if findstr('QA',varnamej)~=0 
            hold on, plot(Tplot(TimeLineQA),bvarplot(is+3:4:end),'r','Linewidth',2) 
        else
            hold on, plot(Tplot(is:end),bvarplot(is:end),'r','Linewidth',2) 
        end
    end
    % GM forecast conditional on exogenous variables paths
    if exopath_yes ==1 
        % Exogenous Paths uncoditional forecast, mean value
        exopplot        = [actual_values; ECfrcst; forecastedexop(:,j,1)];
        % hold on, plot(T(is:end),exopplot(is:end),'g','Linewidth',2) 
        if findstr('QA',varnamej)~=0 
            hold on, plot(Tplot(TimeLineQA),exopplot(is+3:4:end),'r','Linewidth',2) 
        else
            hold on, plot(Tplot(is:end),exopplot(is:end),'r','Linewidth',2) 
        end
    end   
    % GM forecast conditional on endogenous variables paths
    if endopath_yes ==1 
        % Endogenous Paths forecast, mean value
        endopplot        = [actual_values; ECfrcst; forecastedendop(:,j,1)];
        % hold on, plot(T(is:end),endopplot(is:end),'Linewidth',2,'Color',[0 0.5 0]) 
        if findstr('QA',varnamej)~=0 
            hold on, plot(Tplot(TimeLineQA),endopplot(is+3:4:end),'Linewidth',2,'Color',[0 0.5 0]) 
        else
            hold on, plot(Tplot(is:end),endopplot(is:end),'Linewidth',2,'Color',[0 0.5 0]) 
        end
    end   
    % ECFIN forecast
    if frcst_ext == 0
        if ecfin_yes == 1,
            % hold on, plot(T(is:end),EC_plot(is:end),'k-.','Linewidth',2)
            if findstr('QA',varnamej)~=0
                hold on, plot(Tplot(TimeLineQA),EC_plot(is+3:4:end),'k-.','Linewidth',2)               
            else
                hold on, plot(Tplot(is:end),EC_plot(is:end),'k-.','Linewidth',2)
            end
            
            
        end
    end

    if exog_assmpt_yes ==1
        % forecast with assumptions
        exogassmptplot = [smoothvar_exog_assmpt(:,j,1)];
        if findstr('QA',varnamej)~=0 
            hold on, plot(Tplot(end-forecast_:end),exogassmptplot(end-forecast_:end),'g','Linewidth',2) 
        else
            hold on, plot(Tplot(end-forecast_:end),exogassmptplot(end-forecast_:end),'g','Linewidth',2) 
        end
    end
    
    if exog_assmpt_yes ==2
        % forecast with assumptions
        exogassmptplot = [smoothvar_exog_assmpt(:,j,1)];
        if findstr('QA',varnamej)~=0 
            hold on, plot(Tplot(end-forecast_:end),exogassmptplot(end-forecast_:end),'g','Linewidth',2) 
        else
            hold on, plot(Tplot(end-forecast_:end),exogassmptplot(end-forecast_:end),'g','Linewidth',2) 
        end
    end
    
    
    if assmpt_range_yes == 1
              
        for jj=1+(jjj-1)*assmpt_range:jjj*assmpt_range
            assmpt_range_plot = smoothvar_range_assmpt(:,j,jj);
            if assmpt_range == 2
                if rem(jj,2)~=0
                    hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'--+r','Linewidth',1) 
                else
                    hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'--r','Linewidth',1)
                end
            else
                hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'--r','Linewidth',1)
            end
        end
        
    end

    
    if assmpt_range_yes == 2
        
        %  M_.jrc.assmpt_range_vars{jj0} '_' M_.jrc.assmpt_range_vars{jj2}        
        
        
        for jj=1+(jjj-1)*length(M_.jrc.assmpt_range_vars):jjj*length(M_.jrc.assmpt_range_vars)
            %for j0 = 1:length(M_.jrc.assmpt_range_vars)-1
            %for j2 = j0+1:length(M_.jrc.assmpt_range_vars)
            assmpt_range_plot = smoothvar_range_assmpt(:,j,jj);
            if jj == 1+(jjj-1)*length(M_.jrc.assmpt_range_vars)
                %tmp_name = 'HIGH_HIGH';
                hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'-r','Linewidth',2)
            elseif jj == 1+(jjj-1)*length(M_.jrc.assmpt_range_vars)+1
                %tmp_name = 'HIGH_LOW';
                hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'--m','Linewidth',2)
            elseif jj == 1+(jjj-1)*length(M_.jrc.assmpt_range_vars)+2
                %tmp_name = 'LOW_HIGH';
                hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),':c','Linewidth',2)
            elseif jj == 1+(jjj-1)*length(M_.jrc.assmpt_range_vars)+3
                %tmp_name = 'LOW_LOW';
                hold on, plot(Tplot(end-forecast_:end),assmpt_range_plot(end-forecast_:end),'-.k','Linewidth',2)
            end
        end
        
    end
        
        
    
    % GM unconditional forecasts
    % hold on, plot(T(is:end),vplot(is:end),'b','Linewidth',2);
    if findstr('QA',varnamej)~=0 
        if frcst_ext == 1
        %hold on, plot(Tplot(is:4:endQA),vplot1(is+3:4:end),'b','Linewidth',2), plot(Tplot(is:4:endQA),vplot2(is+3:4:end),'k','Linewidth',2), hold on, plot(Tplot(is:4:endQA), vplot3(is+3:4:end),'k--','Linewidth',2)   
        hold on, plot(Tplot(TimeLineQA),vplot1_p,'b','Linewidth',2), plot(Tplot(TimeLineQA),vplot2_p,'k','Linewidth',2), hold on, plot(Tplot(TimeLineQA), vplot3_p,'k--','Linewidth',2)   
        else
        hold on, plot(Tplot(TimeLineQA),vplot(is+3:4:end),'b','Linewidth',2)              
        end
    else
        if frcst_ext == 1
          hold on, plot(Tplot(is:end),vplot1(is:end),'b','Linewidth',2), plot(Tplot(is:end),vplot2(is:end),'k','Linewidth',2), hold on, plot(Tplot(is:end), vplot3(is:end),'k--','Linewidth',2) 
        else
          hold on, plot(Tplot(is:end),vplot(is:end),'b','Linewidth',2)
        end
    end

    
    title(pplottitles{j},'interpreter','none')
    if findstr('QA',varnamej)~=0     
        hold on, plot([Tplot(is) Tplot(end-3)],[steadyvar(j) steadyvar(j)],'k:')
    else
        hold on, plot([Tplot(is) Tplot(end)],[steadyvar(j) steadyvar(j)],'k:')
    end
    axis tight
    
    
    vplot0         = forecastedvar(:,j,1);
    ycomp          = forecastedvar(:,j,1)*NaN;
    yplot          = vplot0;
    yplot_sterr    = [forecastedvar(:,j,2) forecastedvar(:,j,3)];
    ysteady        = repmat(steadyvar(j),[length(TF) 1]);
    xlsinfo(1,1:6) = {'time', pplotvar{j}, '90% conf.', '90% conf.', 'steady state' , ['Comp ',pplotvar{j}]};
    xlsinfo(2:length(TF)+1,1:6) = num2cell([TF yplot yplot_sterr ysteady ycomp]);
   % xlswrite([M_.fname,'_frcst.xls'],xlsinfo,pplotvar{j});
    
    if jplot==nbofplots(jfig) || j==length(pplotvar)
            if assmpt_range_yes == 2
                annotation('textbox', [0.05,0,0.2,0.05],'String', 'HIGH HIGH','Color','r','horizontalalignment','center','interpreter','none');
                annotation('textbox', [0.3,0,0.2,0.05],'String', 'HIGH LOW','Color','m','horizontalalignment','center','interpreter','none');
                annotation('textbox', [0.55,0,0.2,0.05],'String', 'LOW HIGH','Color','c','horizontalalignment','center','interpreter','none');
                annotation('textbox', [0.8,0,0.2,0.05],'String', 'LOW LOW','Color','k','horizontalalignment','center','interpreter','none');                
            end
        files(jfig).name = [M_.fname fnam_suffix '_' int2str(jfig)];
        dyn_saveas(gcf,files(jfig).name,options_.nodisplay,options_.graph_format);
        jplot=0;
    end
end
jplot = 0;
jfig = 0;
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


