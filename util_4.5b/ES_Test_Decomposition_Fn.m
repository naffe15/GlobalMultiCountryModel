function ES_Test_Decomposition_Fn(M_)
global options_

load data T

GPOPpos=strmatch('W_GPOP', M_.param_names);
GP0pos=strmatch('GP0', M_.param_names);
GTFPpos=strmatch('W_GTFP', M_.param_names);
TVATpos=strmatch('ES_TVAT', M_.param_names);
TCpos=strmatch('ES_TC', M_.param_names);
DELTAEpos=strmatch('DELTAE', M_.param_names);
SSCpos=strmatch('SSC', M_.param_names);


W_GPOP=M_.params(GPOPpos);
W_GTFP=M_.params(GTFPpos);
GP0=M_.params(GP0pos);
ES_TVAT=M_.params(TVATpos);
ES_TC=M_.params(TCpos);
DELTAE=M_.params(DELTAEpos);
SSC=M_.params(SSCpos);




wgpop=W_GPOP;
wgtfp=W_GTFP;
gpo=GP0;


% E_TBYN = (exp(E_LPX)/exp(E_LPY)*exp(E_LEX)-exp(E_LPM)/exp(E_LPY)*exp(E_LIM))/exp(E_LY) + E_ZEPS_TB;
% a1=get_smooth('E_LPX');
% c2=-get_smooth('E_LPY');
% c3=get_smooth('E_LEX');
% c4=-get_smooth('E_LPM');
% c5=-get_smooth('E_LIM');
% c6=-get_smooth('E_LY');
% c7=-get_smooth('E_LIM');



% YY0=get_smooth('E_TBYN');

% plot([YY0 c1 c2 c3 c4]),



% % E_LISN=E_LI-E_LY+E_LPI-E_LPY;
% c1=get_smooth('E_LI');
% c2=-get_smooth('E_LY');
% c3=get_smooth('E_LPI');
% c4=-get_smooth('E_LPY');
% YY0=get_smooth('E_LISN');
% plot([YY0 c1 c2 c3 c4]),
% figure, plot([YY0-(c1+c2+c3+c4)])

% E_LGY=E_LG-E_LY;
% c1=get_smooth('E_LG');
% c2=get_smooth('E_LY');
% YY0=get_smooth('E_LGY');
% plot([YY0 c1 c2])


% % E_LIHOUSESN = E_LIHOUSEY + E_LPHOUSEPY;
% c1=get_smooth('E_LIHOUSEY');
% c2=get_smooth('E_LPHOUSEPY');
% YY0=get_smooth('E_LIHOUSESN');
% plot([YY0 c1 c2])

% % E_LCSN = E_LCY + E_LPC - E_LPY;
elpy=get_smooth('ES_LPY');
elcy=get_smooth('ES_LCY');
elpc=get_smooth('ES_LPC');
% c1=elcy;
% c2=elpc;
% c3=elpy;
% YY0=get_smooth('E_LCSN');
% plot([YY0 c1 c2 -c3]), legend({'GBY' 'Others' 'ZEPS_TAX'})



einomgov=get_smooth('ES_INOMGOV');
elb=get_smooth('ES_LBR')+elpy;

elpc=get_smooth('ES_LPC');
elpconstr=get_smooth('ES_LPCONSTR');
elpg=get_smooth('ES_LPG');
elg=get_smooth('ES_LG');
elig=get_smooth('ES_LIG');
% eltr=get_smooth('E_LTR');
eltrobs=get_smooth('ES_LTROBS');
% eben=get_smooth('E_BEN');
ell=get_smooth('ES_LL');
etvat=ES_TVAT;
elc=get_smooth('ES_LC');
elphouse=get_smooth('ES_LPHOUSE');
elihouse=get_smooth('ES_LIHOUSE');
etc=ES_TC;
elk=get_smooth('ES_LK');
%edeltae=get_smooth('ES_DELTAE');
edeltae = DELTAE;
elpi=get_smooth('ES_LPI');
ely=get_smooth('ES_LY');
elwr=get_smooth('ES_LWR');
etl=get_smooth('ES_TL');
% etax=get_smooth('E_TAX')*0;
ezepstax=get_smooth('ES_ZEPS_TAX')-get_smooth('ES_ZEPS_GBY');
ssc=SSC;



C1=    -(einomgov(1:end-1))/(1+wgpop+wgtfp+gpo).*exp(elb(1:end-1))./exp(elpy(2:end));
C1bis=    exp(elb(1:end-1))./exp(elpy(2:end));

C1y=C1./exp(ely(2:end));
C1ybis=C1bis./exp(ely(2:end));
%      -(E_INOM(-1)+W_GPOP+W_GTFP+GP0)*exp(E_LB(-1))/exp(E_LPY)
l1='LB';

C21=   -exp(elpg(2:end))./exp(elpy(2:end)).*exp(elg(2:end));
C21y=C21./exp(ely(2:end));
%      -exp(E_LPC)/exp(E_LPY)*exp(E_LG)
l21='LG';

C22  = -exp(elpconstr(2:end))./exp(elpy(2:end)).*exp(elig(2:end));
C22y=C22./exp(ely(2:end));
%      -exp(E_LPC)/exp(E_LPY)*exp(E_LIG)
l22='LIG';

% C3  = -exp(eltr(2:end))./exp(elpy(2:end))-eben(2:end)./exp(elpy(2:end)).*(1-exp(ell(2:end)));
C3  = -exp(eltrobs(2:end))./exp(elpy(2:end));
C3y=C3./exp(ely(2:end));
%     -exp(E_LTR)/exp(E_LPY)-E_BEN/exp(E_LPY)*(1-exp(E_LL))
l3='LTR and BEN';
C4 =  +etvat.*exp(elpc(2:end))./exp(elpy(2:end)).*exp(elc(2:end))+etvat.*exp(elphouse(2:end))./exp(elpy(2:end)).*exp(elihouse(2:end));
C4y=C4./exp(ely(2:end));
%     +E_TVAT*exp(E_LPC)/exp(E_LPY)*exp(E_LC)+E_TVAT*exp(E_LPHOUSE)/exp(E_LPY)*exp(E_LIHOUSE)
l4='VAT';
C5 =  +etc.*(exp(ely(2:end))-exp(elwr(2:end)).*exp(ell(2:end)))-etc.*edeltae.*exp(elk(2:end)).*exp(elpi(2:end))./exp(elpy(2:end));
C5y=C5./exp(ely(2:end));
%     +E_TC*(exp(E_LY)-exp(E_LWR)*exp(E_LL))-E_TC*E_DELTAE*exp(E_LK)*exp(E_LPI)/exp(E_LPY)
l5='TC';


%    tw

C6 =+(etl(2:end)+ssc).*exp(elwr(2:end)).*exp(ell(2:end));
C6y=C6./exp(ely(2:end));
%   +(E_TL+SSC)*exp(E_LWR)*exp(E_LL)
l6='TL';
% C7=  +(etax(2:end))./exp(elpy(2:end));
% C7y=C7./exp(ely(2:end));
% l7='TAX';
%    +(E_TAX)/exp(E_LPY)
% C8= +(ezepstax(2:end))./exp(elpy(2:end));
% C8y=C8./exp(ely(2:end));
C8= +(ezepstax(2:end));
C8y=C8;
l8='ZEPSTAX';
%   +(E_ZEPS_TAX)/exp(E_LPY)

ss1 = -(get_mean('ES_INOMGOV'))/(1+W_GPOP+W_GTFP+GP0)*get_mean('ES_B');
ss1bis = get_mean('ES_B');
ss21=-get_mean('ES_G');
ss22=-get_mean('ES_IG');
ss3 = -get_mean('ES_TROBS');
ss4=ES_TVAT*get_mean('ES_C')+ES_TVAT*get_mean('ES_IHOUSE');
ss5=ES_TC*(1-get_mean('ES_WS'))-ES_TC*DELTAE*get_mean('ES_K');
ss6=(get_mean('ES_TL')+ssc)*get_mean('ES_WS');
ss8=0;
SS=[ss1 ss3 ss21 ss22 ss6 ss4 ss5 ss8];
SSbis=[ss1bis ss3 ss21 ss22 ss6 ss4 ss5 ss8];

C=[C1y C3y C21y C22y C6y C4y C5y C8y];
Cbis=[C1ybis C3y C21y C22y C6y C4y C5y C8y];

YY=C1y+C21y+C22y+C3y+C4y+C5y+C6y+C8y;
Y=C1+C21+C22+C3+C4+C5+C6+C8;

YY0=get_smooth('ES_GBY');
ll={l1 l3 l21 l22 l6 l4 l5 l8};

% plot(YY,'LineWidth',6)
% hold on, plot(C), legend(ll)

ipos=C>0;
ineg=C<0;

colorsname={'Red'; 'Yellow'; 'Green'; 'Blue'; 'Magenta'; 'Wheat'; 'Lime'; 'Cyan'};
MAP=zeros(length(colorsname),3);

for i=1:length(colorsname)
    MAP(i,:)=rgb(colorsname(i,:));
end

for ii=1:size(C,1)
    TT(ii)=T(options_.first_obs)+0.25+(ii-1)*0.25;

end

figure,
hp=bar((C.*ipos),'stacked');
hold on, hn=bar((C.*ineg),'stacked');
shading faceted; legend(ll); colormap(MAP)
hold on, plot(YY, '-k','LineWidth', 2)
hold on, plot(YY0(2:end), ':k','LineWidth', 2)
set(gca,'Xtick',[4:4:length(C)])
set(gca,'Xticklabel',TT(4:4:end))
saveas(gcf,'Decomposition.fig')

figure,
for j=1:8,
    subplot(3,3,j),
    dump=abs(Cbis(:,j));
    if j==8, dump=(Cbis(:,j)); end
    plot(TT,dump),
    hold on, plot([TT(1),TT(length(C(:,j)))],abs(SSbis([j j])),'r')
    legend(ll{j}),
end
saveas(gcf,'Decomposition_GBY_components.fig')
  



figure,
plot(TT,[YY YY-C(:,8) C(:,8)]), legend({'GBY' 'Others' 'ZEPS_TAX'})
saveas(gcf,'Decomposition_TAX.fig')

