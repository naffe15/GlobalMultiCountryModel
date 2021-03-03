function GDP_pred_err_decomp(obsGDP,obsGDPlevels,obsname,vv)
global options_ bayestopt_ M_ oo_

for i=1:length(obsGDP)
    itemp=strmatch(obsGDP{i},obsname);
    iGDP(i)=itemp(1);
end

err_pred=[];
for t=1:size(vv,1)
% err_pred(:,t)=[(1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PCVAT_IT')*get_mean('C_IT').*((vv(t,iGDP(2)))));  %GC
%                (1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PCVAT_IT')*get_mean('C_IT').*((vv(t,iGDP(3))))); %PHICVAT
%                (1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PI_IT')*get_mean('I_IT').*((vv(t,iGDP(4))))); %GI
%                (1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PI_IT')*get_mean('I_IT').*((vv(t,iGDP(5))))); %PHII
%                (1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PG_IT')*get_mean('G_IT').*((vv(t,iGDP(6))))); %GG
%                (1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PG_IT')*get_mean('G_IT').*((vv(t,iGDP(7))))); %PHIG
%                (1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PIG_IT')*get_mean('IG_IT').*((vv(t,iGDP(8))))); %GIG
%                (1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PIG_IT')*get_mean('IG_IT').*((vv(t,iGDP(9))))); %PHIIG
%                (1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PX_IT')*get_mean('X_IT').*((vv(t,iGDP(10))))); %GXN
%                -(1/(get_mean('YOBS_IT')*get_mean('PYOBS_IT'))).*(get_mean('PMTOT_IT')*get_mean('MTOT_IT').*((vv(t,iGDP(11))))); %GMTOTN
%                -vv(t,iGDP(12));];  %PHIYOBS
           err_pred(:,t)=[(1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{3})*get_mean(obsGDPlevels{2})).*((vv(t,iGDP(2))));  %GC
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{3})*get_mean(obsGDPlevels{2})).*((vv(t,iGDP(3)))-(vv(t,iGDP(14)))); %PHICVAT-PHIYOBS
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{5})*get_mean(obsGDPlevels{4})).*((vv(t,iGDP(4)))); %GI
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{5})*get_mean(obsGDPlevels{4})).*((vv(t,iGDP(5)))-(vv(t,iGDP(14)))); %PHII-PHIYOBS
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{7})*get_mean(obsGDPlevels{6})).*((vv(t,iGDP(6)))); %GG
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{7})*get_mean(obsGDPlevels{6})).*((vv(t,iGDP(7)))-(vv(t,iGDP(14)))); %PHIG-PHIYOBS
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{9})*get_mean(obsGDPlevels{8})).*((vv(t,iGDP(8)))); %GIG
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{9})*get_mean(obsGDPlevels{8})).*((vv(t,iGDP(9)))-(vv(t,iGDP(14)))); %PHIIG-PHIYOBS
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{11})*get_mean(obsGDPlevels{10})).*((vv(t,iGDP(10)))); %GX
                          (1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{11})*get_mean(obsGDPlevels{10})).*((vv(t,iGDP(11)))-(vv(t,iGDP(14)))); %PHIX-PHIYOBS
                         -(1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{13})*get_mean(obsGDPlevels{12})).*((vv(t,iGDP(12)))); %GMTOT
                         -(1/(get_mean(obsGDPlevels{1})*get_mean(obsGDPlevels{14}))).*(get_mean(obsGDPlevels{13})*get_mean(obsGDPlevels{12})).*((vv(t,iGDP(13)))-(vv(t,iGDP(14))));];   %PHIMTOT-PHIYOBS
end

if max(abs(sum(err_pred)-vv(:,iGDP(1))'))>1e-10
    warning('The prediction error of GDP is not equal to the sum of prediction errors of its components!')
end
ipos=err_pred>0;
ineg=err_pred<=0;
nbcmpts=size(err_pred,1);
func = @(x) colorspace('RGB->Lab',x);
MAP = distinguishable_colors(nbcmpts,'w',func);
MAP(end,:) = [0.7 0.7 0.7];
 h = dyn_figure(options_.nodisplay, 'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','landscape','renderermode','auto');
 colormap(MAP);
 set(gcf,'position' ,[50 50 1100 850])  
bar((ipos.*err_pred)','stacked')
hold on
bar((ineg.*err_pred)','stacked')

legend('GC','PHICVAT-PHIYOBS','GI','PHII-PHIYOBS','GG','PHIG-PHIYOBS','GIG','PHIIG-PHIYOBS','GX','PHIX-PHIYOBS','GMTOT','PHIMTOT-PHIYOBS')

plot(vv(:,1),'k','Linewidth',2)

dyn_saveas(gcf,['GDP_pred_err_decomp'],options_.nodisplay,options_.graph_format);
