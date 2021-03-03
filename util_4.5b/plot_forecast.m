function  plot_forecast(pplotvar, M_, options_,oo_,T, Tlim, nplots, oo2,pplottitles )


TT=T(1):0.25:2030;
TF=T(end)+[0.25:0.25:options_.forecast/4]';
TT = [T; TF];
i0=1;
i1=length(TT);
% Tlim=[2005 2020];
if nargin<6 || isempty(Tlim)
    Tlim=[TT(1) TT(end)];
    Tlim=[2005 2020];
end
if nargin<7 || isempty(nplots)
    nplots=[3 2];
    nfigs = ceil(length(pplotvar)/6);
    nplots = repmat(nplots,nfigs,1);
end

if nargin<8 || isempty(oo2),
    comp_flag=0;
else
    comp_flag=1;
end

if nargin<9 || isempty(pplottitles),
    pplottitles=pplotvar;
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
if comp_flag==0,
    oo2=oo_;
end

if comp_flag
    fnam_suffix = '_smooth_and_comp_frcst_';
else
    fnam_suffix = '_smooth_and_frcst_';
end
delete(['*',fnam_suffix,'*.*']);



for j=1:length(pplotvar)
    forecastedvar(:,j)=eval(['oo_.forecast.Mean.',pplotvar{j},';']);
    forecastedvar2(:,j)=eval(['oo2.forecast.Mean.',pplotvar{j},';']);
    smoothedvar(:,j)=eval(['SmoothedVariables.',pplotvar{j},';']);
    steadyvar(j)=get_mean(pplotvar{j});
end

i1=min(i1,size(forecastedvar,1));

jfig=0;
jplot=0;
nplotx = max(nplots);

for j=1:length(pplotvar)
    xlsinfo={};
    
    varnamej=pplotvar{j};
    varpos= strmatch(varnamej, M_.endo_names, 'exact');
    texname{j,:}=M_.endo_names_tex(varpos,:);
    varname{j,:}=varnamej;
    
    
    if jplot==0,
        dyn_figure(options_.nodisplay,'name','Projection plots');
        jfig=jfig+1;
    end
    jplot=jplot+1;
%     subplot(nplots(jfig,1),nplots(jfig,2),jplot)
    subplot(nplots(1),nplots(2),jplot)
    vplot0=forecastedvar(:,j);
    %     vplot0=eval(['ForecastedVariables.',pplotvar{j},';']);
    vplotStd1=eval(['oo_.forecast.HPDsup.',pplotvar{j},';']);
    vplotStd_1=eval(['oo_.forecast.HPDinf.',pplotvar{j},';']);
    if comp_flag,
        ycomp=forecastedvar2(:,j);
        hold on, plot(TF,ycomp,'k--','linewidth',1),
    else
        ycomp=forecastedvar2(:,j)*NaN;
    end
    yplot=vplot0;
    hold on, plot(TF,yplot,'b'),
    yplot_sterr=[vplotStd1 vplotStd_1];
    hold on, plot(TF,yplot_sterr, 'b:')
    title(pplottitles{j},'interpreter','none')
    vplot=smoothedvar(:,j);
    %     vplot=eval(['SmoothedVariables.',pplotvar{j},';']);
    hold on, plot(T(1:length(vplot)),get_smooth(pplotvar{j}),'r.')
    hold on, plot([TT(1) TT(end)],[steadyvar(j) steadyvar(j)],'k:')
    ysteady = repmat(steadyvar(j),[length(TF) 1]);
    set(gca,'xlim',[Tlim(1)-1 Tlim(2)+0.25])
    xlsinfo(1,1:6) = {'time', pplotvar{j}, '90% conf.', '90% conf.', 'steady state' , ['Comp ',pplotvar{j}]};
    xlsinfo(2:length(TF)+1,1:6)=num2cell([TF yplot yplot_sterr ysteady ycomp]);
    if ~ismac
        xlswrite([M_.fname,'_frcst.xls'],xlsinfo,pplotvar{j});
    else
        xlswrite_MACOS([M_.fname,'_frcst.xls'],xlsinfo,pplotvar{j});
    end
    if jplot==nbofplots(jfig) || j==length(pplotvar)
        dyn_saveas(gcf,[M_.fname,fnam_suffix, int2str(Tlim(2)),'_',int2str(jfig)],options_.nodisplay,options_.graph_format);
        jplot=0;
    end
end

if comp_flag
    fidTeX = fopen([M_.fname '_Forecasts_Comp.TeX'],'w');
else
    fidTeX = fopen([M_.fname '_Forecasts.TeX'],'w');
end
files=dir(['*',fnam_suffix,'*.eps']);

h1=1;

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


