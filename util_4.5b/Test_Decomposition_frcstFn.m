function Test_Decomposition_frcstFn_New(ForecastedVariables, M_,forecastflag,addsteady,oo_)

if nargin < 3
    forecastflag=1;
end
if nargin < 4
    addsteady=0;
end

load data T;

GPOPpos=strmatch('W_GPOP', M_.param_names);
GP0pos=strmatch('GP0', M_.param_names);
GTFPpos=strmatch('W_GTFP', M_.param_names);


W_GPOP=M_.params(GPOPpos);
W_GTFP=M_.params(GTFPpos);
GP0=M_.params(GP0pos);

CommonDivisor=1+W_GTFP+W_GPOP+GP0;


IBYComp=ForecastedVariables.E_IBY+get_mean('E_IBY')*(1-addsteady);
TRYComp=ForecastedVariables.E_TRY+get_mean('E_TRY')*(1-addsteady);
TLYComp=ForecastedVariables.E_TLY+get_mean('E_TLY')*(1-addsteady);
VATYComp=ForecastedVariables.E_VATY+get_mean('E_VATY')*(1-addsteady);
TCYComp=ForecastedVariables.E_TCY+get_mean('E_TCY')*(1-addsteady);
TAXYComp=ForecastedVariables.E_TAXY+get_mean('E_TAXY')*(1-addsteady)-ForecastedVariables.E_ZEPS_GBY;
GSComp=ForecastedVariables.E_GS+get_mean('E_GS')*(1-addsteady);
IGSComp=ForecastedVariables.E_IGS+get_mean('E_IGS')*(1-addsteady);
if forecastflag
    TT=1:length(IGSComp);
else
    TT=1:length(oo_.SmoothedVariables.E_GBY);
end

legendstr={'IBY','TRY','GS','IGS','TLY','VAT','TCY','TAX'};
Components=[-IBYComp -TRYComp -GSComp -IGSComp TLYComp VATYComp TCYComp TAXYComp];
% Components=Components/CommonDivisor;
Components=Components(1:TT(end),:);
GBYSum=sum(Components,2);
%figure('units','normalized','outerposition',[0 0 1 1])
plot([GBYSum oo_.SmoothedVariables.E_GBY+get_mean('E_GBY') ForecastedVariables.E_GBY(1:TT(end))])
legend({'Sum','Smooth' 'Fore'})
% hold all, plot(oo_.SmoothedVariables.E_GBY)
ipos=Components>0;
ineg=Components<0;

colorsname={'Red'; 'Yellow'; 'Green'; 'Blue'; 'Magenta'; 'Wheat'; 'Lime'; 'Cyan'};
MAP=zeros(length(colorsname),3);

for i=1:length(colorsname)
    MAP(i,:)=rgb(colorsname(i,:));
end

for ii=1:size(GBYSum,1)
TT0(ii)=T(1)+0.25+(ii-1)*0.25;

end

h0=figure;
% fullscreen = get(0,'ScreenSize');
% figure('Position',[0 0 fullscreen(3) fullscreen(4)])
% figure('units','normalized','outerposition',[0 0 1 1])
hp=bar((Components.*ipos),'stacked');
hold on, hn=bar((Components.*ineg),'stacked');
shading faceted; colormap(MAP)
hold on, plot([ForecastedVariables.E_GBY(1:TT(end)) GBYSum],'-k','LineWidth', 2)
% hold on, plot(netgov,'-k','LineWidth', 2)
% title(['TGOVB1=',num2str(TGOVB1E),'; TGOVB2=',num2str(TGOVB2E)]),

set(gca,'Xtick',[5:4:length(TT0)])
set(gca,'Xticklabel',TT0(4:4:end))

legend(legendstr,'Location','SouthEast','Orientation','horizontal')
saveas(gcf,'Decomposition_New.fig'),
% plot([ForecastedVariables.E_GBY GBYSum])

