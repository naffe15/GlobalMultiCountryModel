function DE_Test_Decomposition_Fn(M_)
global options_

load data T

GPOPpos=strmatch('W_GPOP', M_.param_names);
GP0pos=strmatch('GP0', M_.param_names);
GTFPpos=strmatch('W_GTFP', M_.param_names);
TVATpos=strmatch('DE_TVAT', M_.param_names);
TCpos=strmatch('DE_TC', M_.param_names);
DELTAEpos=strmatch('DELTAE', M_.param_names);
SSCpos=strmatch('SSC', M_.param_names);


W_GPOP=M_.params(GPOPpos);
W_GTFP=M_.params(GTFPpos);
GP0=M_.params(GP0pos);
DE_TVAT=M_.params(TVATpos);
DE_TC=M_.params(TCpos);
DELTAE=M_.params(DELTAEpos);
SSC=M_.params(SSCpos);




wgpop=W_GPOP;
wgtfp=W_GTFP;
gpo=GP0;


% DE_TBYN = (exp(DE_LPX)/exp(DE_LPY)*exp(DE_LEX)-exp(DE_LPM)/exp(DE_LPY)*exp(DE_LIM))/exp(DE_LY) + DE_ZEPS_TB;
% a1=get_smooth('DE_LPX');
% c2=-get_smooth('DE_LPY');
% c3=get_smooth('DE_LEX');
% c4=-get_smooth('DE_LPM');
% c5=-get_smooth('DE_LIM');
% c6=-get_smooth('DE_LY');
% c7=-get_smooth('DE_LIM');



% YY0=get_smooth('DE_TBYN');

% plot([YY0 c1 c2 c3 c4]),



% % DE_LISN=DE_LI-DE_LY+DE_LPI-DE_LPY;
% c1=get_smooth('DE_LI');
% c2=-get_smooth('DE_LY');
% c3=get_smooth('DE_LPI');
% c4=-get_smooth('DE_LPY');
% YY0=get_smooth('DE_LISN');
% plot([YY0 c1 c2 c3 c4]),
% figure, plot([YY0-(c1+c2+c3+c4)])

% DE_LGY=DE_LG-DE_LY;
% c1=get_smooth('DE_LG');
% c2=get_smooth('DE_LY');
% YY0=get_smooth('DE_LGY');
% plot([YY0 c1 c2])


% % DE_LIHOUSESN = DE_LIHOUSEY + DE_LPHOUSEPY;
% c1=get_smooth('DE_LIHOUSEY');
% c2=get_smooth('DE_LPHOUSEPY');
% YY0=get_smooth('DE_LIHOUSESN');
% plot([YY0 c1 c2])

% % DE_LCSN = DE_LCY + DE_LPC - DE_LPY;
elpy=get_smooth('DE_LPY');
elcy=get_smooth('DE_LCY');
elpc=get_smooth('DE_LPC');
% c1=elcy;
% c2=elpc;
% c3=elpy;
% YY0=get_smooth('DE_LCSN');
% plot([YY0 c1 c2 -c3]), legend({'GBY' 'Others' 'ZEPS_TAX'})



einomgov=get_smooth('DE_INOMGOV');
elb=get_smooth('DE_LBR')+elpy;

elpc=get_smooth('DE_LPC');
elpconstr=get_smooth('DE_LPCONSTR');
elpg=get_smooth('DE_LPG');
elg=get_smooth('DE_LG');
elig=get_smooth('DE_LIG');
% eltr=get_smooth('DE_LTR');
eltrobs=get_smooth('DE_LTROBS');
% eben=get_smooth('DE_BEN');
ell=get_smooth('DE_LL');
etvat=DE_TVAT;
elc=get_smooth('DE_LC');
elphouse=get_smooth('DE_LPHOUSE');
elihouse=get_smooth('DE_LIHOUSE');
etc=DE_TC;
elk=get_smooth('DE_LK');
elpi=get_smooth('DE_LPI');
ely=get_smooth('DE_LY');
elwr=get_smooth('DE_LWR');
etl=get_smooth('DE_TL');
% etax=get_smooth('DE_TAX')*0;
ezepstax=get_smooth('DE_ZEPS_TAX')-get_smooth('DE_ZEPS_GBY');
ssc=SSC;



C1=    -(einomgov(1:end-1))/(1+wgpop+wgtfp+gpo).*exp(elb(1:end-1))./exp(elpy(2:end));
C1bis=    exp(elb(1:end-1))./exp(elpy(2:end));

C1y=C1./exp(ely(2:end));
C1ybis=C1bis./exp(ely(2:end));
%      -(DE_INOM(-1)+W_GPOP+W_GTFP+GP0)*exp(DE_LB(-1))/exp(DE_LPY)
l1='LB';

C21=   -exp(elpg(2:end))./exp(elpy(2:end)).*exp(elg(2:end));
C21y=C21./exp(ely(2:end));
%      -exp(DE_LPC)/exp(DE_LPY)*exp(DE_LG)
l21='LG';

C22  = -exp(elpconstr(2:end))./exp(elpy(2:end)).*exp(elig(2:end));
C22y=C22./exp(ely(2:end));
%      -exp(DE_LPC)/exp(DE_LPY)*exp(DE_LIG)
l22='LIG';

% C3  = -exp(eltr(2:end))./exp(elpy(2:end))-eben(2:end)./exp(elpy(2:end)).*(1-exp(ell(2:end)));
C3  = -exp(eltrobs(2:end))./exp(elpy(2:end));
C3y=C3./exp(ely(2:end));
%     -exp(DE_LTR)/exp(DE_LPY)-DE_BEN/exp(DE_LPY)*(1-exp(DE_LL))
l3='LTR and BEN';
C4 =  +etvat.*exp(elpc(2:end))./exp(elpy(2:end)).*exp(elc(2:end))+etvat.*exp(elphouse(2:end))./exp(elpy(2:end)).*exp(elihouse(2:end));
C4y=C4./exp(ely(2:end));
%     +DE_TVAT*exp(DE_LPC)/exp(DE_LPY)*exp(DE_LC)+DE_TVAT*exp(DE_LPHOUSE)/exp(DE_LPY)*exp(DE_LIHOUSE)
l4='VAT';
C5 =  +etc.*(exp(ely(2:end))-exp(elwr(2:end)).*exp(ell(2:end)))-etc.*DELTAE.*exp(elk(2:end)).*exp(elpi(2:end))./exp(elpy(2:end));
C5y=C5./exp(ely(2:end));
%     +DE_TC*(exp(DE_LY)-exp(DE_LWR)*exp(DE_LL))-DE_TC*DE_DELTAE*exp(DE_LK)*exp(DE_LPI)/exp(DE_LPY)
l5='TC';


%    tw

C6 =+(etl(2:end)+ssc).*exp(elwr(2:end)).*exp(ell(2:end));
C6y=C6./exp(ely(2:end));
%   +(DE_TL+SSC)*exp(DE_LWR)*exp(DE_LL)
l6='TL';
% C7=  +(etax(2:end))./exp(elpy(2:end));
% C7y=C7./exp(ely(2:end));
% l7='TAX';
%    +(DE_TAX)/exp(DE_LPY)
% C8= +(ezepstax(2:end))./exp(elpy(2:end));
% C8y=C8./exp(ely(2:end));
C8= +(ezepstax(2:end));
C8y=C8;
l8='ZEPSTAX';
%   +(DE_ZEPS_TAX)/exp(DE_LPY)

ss1 = -(get_mean('DE_INOMGOV'))/(1+W_GPOP+W_GTFP+GP0)*get_mean('DE_B');
ss1bis = get_mean('DE_B');
ss21=-get_mean('DE_G');
ss22=-get_mean('DE_IG');
ss3 = -get_mean('DE_TROBS');
ss4=DE_TVAT*get_mean('DE_C')+DE_TVAT*get_mean('DE_IHOUSE');
ss5=DE_TC*(1-get_mean('DE_WS'))-DE_TC*DELTAE*get_mean('DE_K');
ss6=(get_mean('DE_TL')+ssc)*get_mean('DE_WS');
ss8=0;
SS=[ss1 ss3 ss21 ss22 ss6 ss4 ss5 ss8];
SSbis=[ss1bis ss3 ss21 ss22 ss6 ss4 ss5 ss8];

C=[C1y C3y C21y C22y C6y C4y C5y C8y];
Cbis=[C1ybis C3y C21y C22y C6y C4y C5y C8y];

YY=C1y+C21y+C22y+C3y+C4y+C5y+C6y+C8y;
Y=C1+C21+C22+C3+C4+C5+C6+C8;

YY0=get_smooth('DE_GBY');
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

