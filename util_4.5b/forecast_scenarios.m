function forecast_scenarios(M_,options_,oo_,T,nfrcs,outx,out,TLim)

sname=outx.sname;
tname=outx.tname;
xname=outx.xname;
txname=outx.txname;
xvalues=outx.xvalues;
dyn_figure(options_.nodisplay,'name','Forecast Scenarios');

TT=T(1):0.25:(T(1)+(length(T)+nfrcs-1)*0.25);

xlsinfo={};
xlsinfo(1,1)={'time'};
xlsinfo(2:length(TT)+1,1)=num2cell(TT');
for j=1:length(sname),
    xlstmp={};
    if length(sname)>6,
        subplot(3,3,j),
    else
        subplot(3,2,j),
    end
    eval(['plot(TT,[outx.ForecastedVariables.',sname{j},'],''k'')']);
    hold all,
    eval(['plot(TT,[out.ForecastedVariables.',sname{j},'],'':k'')']);
    xlstmp(1,1:2)={sname{j},[sname{j},'_base']};
    xlstmp(2:length(TT)+1,1:2)=num2cell([getfield(out.ForecastedVariables,sname{j}) getfield(outx.ForecastedVariables,sname{j})]);
    xlsinfo = [xlsinfo xlstmp];
    plot(T,get_smooth(sname{j}),'.r')
    title(tname{j})
    set(gca,'xlim',TLim)
end
xlswrite([M_.fname,'_frcst.xls'],xlsinfo,'scenarios');
dyn_saveas(gcf,[M_.fname '_frcst_scenarios'],options_.nodisplay,options_.graph_format);

dyn_figure(options_.nodisplay,'name','Forecast Scenarios Shocks');

xlsinfo={};
xlsinfo(1,1)={'time'};
xlsinfo(2:length(TT)+1,1)=num2cell(TT');
for j=1:length(xname),
    xlstmp={};
    if length(xname)>6,
        subplot(3,3,j),
    else
        subplot(3,2,j),
    end
    yval=[getfield(oo_.SmoothedShocks,xname{j}); xvalues(:,j)];
    plot(TT,yval,'k');
    hold all,
    yval0=[getfield(oo_.SmoothedShocks,xname{j}); xvalues(:,j)*0];
    plot(TT,yval0,':k');
    plot(T,getfield(oo_.SmoothedShocks,xname{j}),'.r')
    title(txname{j})
    set(gca,'xlim',TLim)
    xlstmp(1,1:2)={xname{j},[xname{j},'_base']};
    xlstmp(2:length(TT)+1,1:2)=num2cell([yval yval0]);
    xlsinfo = [xlsinfo xlstmp];
end
xlswrite([M_.fname,'_frcst.xls'],xlsinfo,'scenarios_shocks');
dyn_saveas(gcf,[M_.fname '_frcst_scenarios_shocks'],options_.nodisplay,options_.graph_format);


