function  Plot_Forecasts_Fn(pplotvar, M_, options_,oo_, ForecastedVariables,ForecastedVariablesStd,pplotvar0,T, Tlim, nplots, ForecastedVariables2,pplottitles )


TT=T(1):0.25:2030;
i0=1;
i1=length(TT);
% Tlim=[2005 2020];
if nargin<9 || isempty(Tlim)
    Tlim=[TT(1) TT(end)];
    Tlim=[2005 2020];
end
if nargin<10 || isempty(nplots)
    nplots=[3 2];
    nfigs = ceil(length(pplotvar)/6);
    nplots = repmat(nplots,nfigs,1);
end

if nargin<11 || isempty(ForecastedVariables2),
    comp_flag=0;
else
    comp_flag=1;
end

if nargin<12 || isempty(pplottitles),
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
    ForecastedVariables2=ForecastedVariables;
end

if comp_flag
    fnam_suffix = '_smooth_and_comp_frcst_';
else
    fnam_suffix = '_smooth_and_frcst_';
end
delete(['*',fnam_suffix,'*.*']);



for j=1:length(pplotvar)
    forecastedvar(:,j)=eval(['ForecastedVariables.',pplotvar{j},';']);
    forecastedvar2(:,j)=eval(['ForecastedVariables2.',pplotvar{j},';']);
    smoothedvar(:,j)=eval(['SmoothedVariables.',pplotvar{j},';']);
    steadyvar(j)=get_mean(pplotvar{j});
    if forecastedvar(1,j)~=smoothedvar(1,j),
        forecastedvar(:,j) = forecastedvar(:,j) - steadyvar(j);
    end
    if forecastedvar2(1,j)~=smoothedvar(1,j),
        forecastedvar2(:,j) = forecastedvar2(:,j) - steadyvar(j);
    end
    for k=1:size(pplotvar0,1)
        if strcmp(pplotvar{j}, pplotvar0{k,1}),
            steadyvar(j) = steadyvar(j) + pplotvar0{k,2};
            forecastedvar(:,j) = forecastedvar(:,j) + pplotvar0{k,2};
            forecastedvar2(:,j) = forecastedvar2(:,j) + pplotvar0{k,2};
            smoothedvar(:,j) = smoothedvar(:,j) + pplotvar0{k,2};
        end
    end
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
    subplot(nplotx(1),nplotx(2),jplot)
    vplot0=forecastedvar(:,j);
    %     vplot0=eval(['ForecastedVariables.',pplotvar{j},';']);
    vplotStd1=eval(['vplot0+2*ForecastedVariablesStd.',pplotvar{j},';']);
    vplotStd_1=eval(['vplot0-2*ForecastedVariablesStd.',pplotvar{j},';']);
    if comp_flag,
        vplot2=forecastedvar2(:,j);
        ycomp=vplot2(i0:i1)+get_mean(pplotvar{j});
        hold on, plot(TT,ycomp,'k--','linewidth',1),
    else
        ycomp=vplot0(i0:i1)*NaN;
    end
    yplot=vplot0(i0:i1)+get_mean(pplotvar{j});
    hold on, plot(TT,yplot,'b'),
    yplot_sterr=[vplotStd1(i0:i1) vplotStd_1(i0:i1)]+get_mean(pplotvar{j});
    hold on, plot(TT,yplot_sterr, 'b:')
    title(pplottitles{j},'interpreter','none')
    vplot=smoothedvar(:,j);
    %     vplot=eval(['SmoothedVariables.',pplotvar{j},';']);
    hold on, plot(T(1:length(vplot)),vplot+get_mean(pplotvar{j}),'r.')
    hold on, plot([TT(1) TT(end)],[steadyvar(j) steadyvar(j)],'k:')
    ysteady = repmat(steadyvar(j),[length(TT) 1]);
    set(gca,'xlim',[Tlim(1)-1 Tlim(2)+0.25])
    xlsinfo(1,1:6) = {'time', pplotvar{j}, '90% conf.', '90% conf.', 'steady state' , ['Comp ',pplotvar{j}]};
    xlsinfo(2:length(TT)+1,1:6)=num2cell([TT' yplot yplot_sterr ysteady ycomp]);
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


