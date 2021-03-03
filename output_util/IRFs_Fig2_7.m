%% load and save IRFs from oo_
% clear all
cd(workdirectory)
wd=pwd;
mkdir([wd '\figures'])

countries={'DE', 'FR', 'IT', 'ES'};

b = countries;
for i= 1:length(countries)
 results.(b{i})=load([countries{i} '\gemc_results.mat']);
end


%%
% ****************************************
% TFP
% ****************************************
T=1:1:40;
figure(1),
subplot(3,3,1)
plot( results.DE.oo_.irfs.LYOBS_DE_EPS_LAYTREND_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LYOBS_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LYOBS_IT_EPS_LAYTREND_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LYOBS_ES_EPS_LAYTREND_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP','FontSize',20)

subplot(3,3,2)
plot( results.DE.oo_.irfs.LC_DE_EPS_LAYTREND_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LC_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LC_IT_EPS_LAYTREND_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LC_ES_EPS_LAYTREND_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Consumption','FontSize',20)

subplot(3,3,3)
plot( results.DE.oo_.irfs.LI_DE_EPS_LAYTREND_DE(1:40)'.*100,'--','LineWidth',3),hold on, plot(results.FR.oo_.irfs.LI_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LI_IT_EPS_LAYTREND_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LI_ES_EPS_LAYTREND_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Investment','FontSize',20)

subplot(3,3,4)
plot(results.DE.oo_.irfs.LN_DE_EPS_LAYTREND_DE(1:40)'.*100,'--','LineWidth',3),hold on, plot(results.FR.oo_.irfs.LN_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LN_IT_EPS_LAYTREND_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LN_ES_EPS_LAYTREND_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Hours','FontSize',20)

subplot(3,3,5)
plot( results.DE.oo_.irfs.LWR_DE_EPS_LAYTREND_DE(1:40)'.*100,'--','LineWidth',3),hold on, plot(results.FR.oo_.irfs.LWR_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LWR_IT_EPS_LAYTREND_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LWR_ES_EPS_LAYTREND_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real wages','FontSize',20)

subplot(3,3,6)
plot( results.DE.oo_.irfs.R_DE_EPS_LAYTREND_DE(1:40)'.*100,'--','LineWidth',3),hold on, plot(results.FR.oo_.irfs.R_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.R_IT_EPS_LAYTREND_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.R_ES_EPS_LAYTREND_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real interest rate','FontSize',20)

subplot(3,3,7)
plot(results.DE.oo_.irfs.PHIYOBS_DE_EPS_LAYTREND_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.PHIYOBS_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.PHIYOBS_IT_EPS_LAYTREND_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.PHIYOBS_ES_EPS_LAYTREND_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP Inflation','FontSize',20)

subplot(3,3,8)
plot(results.DE.oo_.irfs.LRER_DE_EPS_LAYTREND_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LRER_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LRER_IT_EPS_LAYTREND_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LRER_ES_EPS_LAYTREND_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real exchange rate','FontSize',20)

subplot(3,3,9)
h=plot(T,results.DE.oo_.irfs.TBY_DE_EPS_LAYTREND_DE(1:40)'.*100,'--', T, results.FR.oo_.irfs.TBY_FR_EPS_LAYTREND_FR(1:40)'.*100,'-.', T, results.IT.oo_.irfs.TBY_IT_EPS_LAYTREND_IT(1:40)'.*100,  T, results.ES.oo_.irfs.TBY_ES_EPS_LAYTREND_ES(1:40)'.*100,':', T,0*(1:40),':k','LineWidth',3);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
set(h(5),'linewidth',1);
title('Trade balance to GDP','FontSize',20)
legend(h([1 2 3 4]),{'DE','FR','IT','ES'},'FontSize',20,'Orientation','horizontal','Position',[0.265 0.01 0.5 0.05])
set(gcf, 'Position', get(0, 'Screensize'));
savefig([wd '\figures\' 'TFP'])
print('-f1', '-depsc', [wd '\figures\' 'TFP.eps'])
close(figure(1));

%%
% ****************************************
% Gov. consumption
% ****************************************
figure(2),
subplot(3,3,1)
plot(results.DE.oo_.irfs.LYOBS_DE_EPS_G_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LYOBS_FR_EPS_G_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LYOBS_IT_EPS_G_IT(1:40)'.*100, 'LineWidth', 3), hold on, plot(results.ES.oo_.irfs.LYOBS_ES_EPS_G_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP','FontSize',20)

subplot(3,3,2)
plot(results.DE.oo_.irfs.LC_DE_EPS_G_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LC_FR_EPS_G_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LC_IT_EPS_G_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LC_ES_EPS_G_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Consumption','FontSize',20)

subplot(3,3,3)
plot(results.DE.oo_.irfs.LI_DE_EPS_G_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LI_FR_EPS_G_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LI_IT_EPS_G_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LI_ES_EPS_G_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Investment','FontSize',20)

subplot(3,3,4)
plot(results.DE.oo_.irfs.LN_DE_EPS_G_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LN_FR_EPS_G_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LN_IT_EPS_G_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LN_ES_EPS_G_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Hours','FontSize',20)

subplot(3,3,5)
plot(results.DE.oo_.irfs.LWR_DE_EPS_G_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LWR_FR_EPS_G_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LWR_IT_EPS_G_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LWR_ES_EPS_G_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real wages','FontSize',20)

subplot(3,3,6)
plot(results.DE.oo_.irfs.R_DE_EPS_G_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.R_FR_EPS_G_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.R_IT_EPS_G_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.R_ES_EPS_G_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real interest rate','FontSize',20)

subplot(3,3,7)
plot(results.DE.oo_.irfs.PHIYOBS_DE_EPS_G_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.PHIYOBS_FR_EPS_G_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.PHIYOBS_IT_EPS_G_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.PHIYOBS_ES_EPS_G_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP Inflation','FontSize',20)

subplot(3,3,8)
plot(results.DE.oo_.irfs.LRER_DE_EPS_G_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LRER_FR_EPS_G_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LRER_IT_EPS_G_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LRER_ES_EPS_G_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real exchange rate','FontSize',20)

subplot(3,3,9)
h=plot(T,results.DE.oo_.irfs.TBY_DE_EPS_G_DE(1:40)'.*100,'--', T, results.FR.oo_.irfs.TBY_FR_EPS_G_FR(1:40)'.*100,'-.', T, results.IT.oo_.irfs.TBY_IT_EPS_G_IT(1:40)'.*100,  T, results.ES.oo_.irfs.TBY_ES_EPS_G_ES(1:40)'.*100,':', T,0*(1:40),':k','LineWidth',3);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
set(h(5),'linewidth',1);
title('Trade balance to GDP','FontSize',20)
legend(h([1 2 3 4]),{'DE','FR','IT','ES'},'FontSize',20,'Orientation','horizontal','Position',[0.265 0.01 0.5 0.05])
set(gcf, 'Position', get(0, 'Screensize'));
savefig([wd '\figures\' 'Gov_exp'])
print('-f2', '-depsc', [wd '\figures\' 'Gov_exp.eps'])
close(figure(2));

%%
% ****************************************
% Monetary
% ****************************************
figure(3),
subplot(3,3,1)
plot(results.DE.oo_.irfs.LYOBS_DE_EPS_INOM_EA(1:40)'.*-25,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LYOBS_FR_EPS_INOM_EA(1:40)'.*-25,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LYOBS_IT_EPS_INOM_EA(1:40)'.*-25, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LYOBS_ES_EPS_INOM_EA(1:40)'.*-25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP','FontSize',20)

subplot(3,3,2)
plot(results.DE.oo_.irfs.LC_DE_EPS_INOM_EA(1:40)'.*-25,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LC_FR_EPS_INOM_EA(1:40)'.*-25,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LC_IT_EPS_INOM_EA(1:40)'.*-25, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LC_ES_EPS_INOM_EA(1:40)'.*-25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Consumption','FontSize',20)

subplot(3,3,3)
plot(results.DE.oo_.irfs.LI_DE_EPS_INOM_EA(1:40)'.*-25,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LI_FR_EPS_INOM_EA(1:40)'.*-25,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LI_IT_EPS_INOM_EA(1:40)'.*-25, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LI_ES_EPS_INOM_EA(1:40)'.*-25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Investment','FontSize',20)

subplot(3,3,4)
plot(results.DE.oo_.irfs.LN_DE_EPS_INOM_EA(1:40)'.*-25,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LN_FR_EPS_INOM_EA(1:40)'.*-25,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LN_IT_EPS_INOM_EA(1:40)'.*-25, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LN_ES_EPS_INOM_EA(1:40)'.*-25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Hours','FontSize',20)

subplot(3,3,5)
plot(results.DE.oo_.irfs.LWR_DE_EPS_INOM_EA(1:40)'.*-25,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LWR_FR_EPS_INOM_EA(1:40)'.*-25,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LWR_IT_EPS_INOM_EA(1:40)'.*-25, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LWR_ES_EPS_INOM_EA(1:40)'.*-25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real wages','FontSize',20)

subplot(3,3,6)
plot(results.DE.oo_.irfs.R_DE_EPS_INOM_EA(1:40)'.*-25,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.R_FR_EPS_INOM_EA(1:40)'.*-25,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.R_IT_EPS_INOM_EA(1:40)'.*-25, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.R_ES_EPS_INOM_EA(1:40)'.*-25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real interest rate','FontSize',20)

subplot(3,3,7)
plot(results.DE.oo_.irfs.PHIYOBS_DE_EPS_INOM_EA(1:40)'.*-25,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.PHIYOBS_FR_EPS_INOM_EA(1:40)'.*-25,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.PHIYOBS_IT_EPS_INOM_EA(1:40)'.*-25, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.PHIYOBS_ES_EPS_INOM_EA(1:40)'.*-25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP Inflation','FontSize',20)

subplot(3,3,8)
plot(results.DE.oo_.irfs.LRER_DE_EPS_INOM_EA(1:40)'.*-25,'--','LineWidth',3), hold on, plot(results.FR.oo_.irfs.LRER_FR_EPS_INOM_EA(1:40)'.*-25,'-.','LineWidth',3), hold on, plot(results.IT.oo_.irfs.LRER_IT_EPS_INOM_EA(1:40)'.*-25, 'LineWidth', 2), hold on, plot(results.ES.oo_.irfs.LRER_ES_EPS_INOM_EA(1:40)'.*-25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real exchange rate','FontSize',20)

subplot(3,3,9)
h=plot(T,results.DE.oo_.irfs.TBY_DE_EPS_INOM_EA(1:40)'.*-25,'--', T, results.FR.oo_.irfs.TBY_FR_EPS_INOM_EA(1:40)'.*-25,'-.', T, results.IT.oo_.irfs.TBY_IT_EPS_INOM_EA(1:40)'.*-25,  T, results.ES.oo_.irfs.TBY_ES_EPS_INOM_EA(1:40)'.*-25,':', T,0*(1:40),':k','LineWidth',3);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
set(h(5),'linewidth',1);
title('Trade balance to GDP','FontSize',20)
legend(h([1 2 3 4]),{'DE','FR','IT','ES'},'FontSize',20,'Orientation','horizontal','Position',[0.265 0.01 0.5 0.05])
set(gcf, 'Position', get(0, 'Screensize'));
savefig([wd '\figures\' 'Monetary'])
print('-f3', '-depsc', [wd '\figures\' 'Monetary.eps'])
close(figure(3));

%%
% ****************************************
% Saving shock
% ****************************************

figure(5),
subplot(3,3,1)
plot(-results.DE.oo_.irfs.LYOBS_DE_EPS_UC_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LYOBS_FR_EPS_UC_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LYOBS_IT_EPS_UC_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LYOBS_ES_EPS_UC_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP','FontSize',20)

subplot(3,3,2)
plot(-results.DE.oo_.irfs.LC_DE_EPS_UC_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LC_FR_EPS_UC_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LC_IT_EPS_UC_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LC_ES_EPS_UC_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Consumption','FontSize',20)

subplot(3,3,3)
plot(-results.DE.oo_.irfs.LI_DE_EPS_UC_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LI_FR_EPS_UC_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LI_IT_EPS_UC_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LI_ES_EPS_UC_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Investment','FontSize',20)

subplot(3,3,4)
plot(-results.DE.oo_.irfs.LN_DE_EPS_UC_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LN_FR_EPS_UC_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LN_IT_EPS_UC_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LN_ES_EPS_UC_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Hours','FontSize',20)

subplot(3,3,5)
plot(-results.DE.oo_.irfs.LWR_DE_EPS_UC_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LWR_FR_EPS_UC_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LWR_IT_EPS_UC_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LWR_ES_EPS_UC_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real wages','FontSize',20)

subplot(3,3,6)
plot(-results.DE.oo_.irfs.R_DE_EPS_UC_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.R_FR_EPS_UC_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.R_IT_EPS_UC_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.R_ES_EPS_UC_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real interest rate','FontSize',20)

subplot(3,3,7)
plot(-results.DE.oo_.irfs.PHIYOBS_DE_EPS_UC_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.PHIYOBS_FR_EPS_UC_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.PHIYOBS_IT_EPS_UC_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.PHIYOBS_ES_EPS_UC_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP Inflation','FontSize',20)

subplot(3,3,8)
plot(-results.DE.oo_.irfs.LRER_DE_EPS_UC_DE(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LRER_FR_EPS_UC_FR(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LRER_IT_EPS_UC_IT(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LRER_ES_EPS_UC_ES(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real exchange rate','FontSize',20)

subplot(3,3,9)
h=plot(T,-results.DE.oo_.irfs.TBY_DE_EPS_UC_DE(1:40)'.*100,'--', T, -results.FR.oo_.irfs.TBY_FR_EPS_UC_FR(1:40)'.*100,'-.', T, -results.IT.oo_.irfs.TBY_IT_EPS_UC_IT(1:40)'.*100,  T, -results.ES.oo_.irfs.TBY_ES_EPS_UC_ES(1:40)'.*100,':', T,0*(1:40),':k','LineWidth',3);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
set(h(5),'linewidth',1);
title('Trade balance to GDP','FontSize',20)
legend(h([1 2 3 4]),{'DE','FR','IT','ES'},'FontSize',20,'Orientation','horizontal','Position',[0.265 0.01 0.5 0.05])
set(gcf, 'Position', get(0, 'Screensize'));
savefig([wd '\figures\' 'Saving'])
print('-f5', '-depsc', [wd '\figures\' 'Saving.eps'])
close(figure(5));

%%
% ****************************************
% Saving shock in RoW
% ****************************************

figure(6),
subplot(3,3,1)
plot(-results.DE.oo_.irfs.LYOBS_DE_EPS_UC_RoW(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LYOBS_FR_EPS_UC_RoW(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LYOBS_IT_EPS_UC_RoW(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LYOBS_ES_EPS_UC_RoW(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
% legend('DE','FR','IT','ES')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP','FontSize',20)

subplot(3,3,2)
plot(-results.DE.oo_.irfs.LC_DE_EPS_UC_RoW(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LC_FR_EPS_UC_RoW(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LC_IT_EPS_UC_RoW(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LC_ES_EPS_UC_RoW(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Consumption','FontSize',20)

subplot(3,3,3)
plot(-results.DE.oo_.irfs.LI_DE_EPS_UC_RoW(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LI_FR_EPS_UC_RoW(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LI_IT_EPS_UC_RoW(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LI_ES_EPS_UC_RoW(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Investment','FontSize',20)

subplot(3,3,4)
plot(-results.DE.oo_.irfs.LN_DE_EPS_UC_RoW(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LN_FR_EPS_UC_RoW(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LN_IT_EPS_UC_RoW(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LN_ES_EPS_UC_RoW(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Hours','FontSize',20)

subplot(3,3,5)
plot(-results.DE.oo_.irfs.LWR_DE_EPS_UC_RoW(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LWR_FR_EPS_UC_RoW(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LWR_IT_EPS_UC_RoW(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LWR_ES_EPS_UC_RoW(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real wages','FontSize',20)

subplot(3,3,6)
plot(-results.DE.oo_.irfs.R_DE_EPS_UC_RoW(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.R_FR_EPS_UC_RoW(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.R_IT_EPS_UC_RoW(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.R_ES_EPS_UC_RoW(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real interest rate','FontSize',20)

subplot(3,3,7)
plot(-results.DE.oo_.irfs.PHIYOBS_DE_EPS_UC_RoW(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.PHIYOBS_FR_EPS_UC_RoW(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.PHIYOBS_IT_EPS_UC_RoW(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.PHIYOBS_ES_EPS_UC_RoW(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP Inflation','FontSize',20)

subplot(3,3,8)
plot(-results.DE.oo_.irfs.LRER_DE_EPS_UC_RoW(1:40)'.*100,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LRER_FR_EPS_UC_RoW(1:40)'.*100,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LRER_IT_EPS_UC_RoW(1:40)'.*100, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LRER_ES_EPS_UC_RoW(1:40)'.*100,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real exchange rate','FontSize',20)

subplot(3,3,9)
h=plot(T,-results.DE.oo_.irfs.TBY_DE_EPS_UC_RoW(1:40)'.*100,'--', T, -results.FR.oo_.irfs.TBY_FR_EPS_UC_RoW(1:40)'.*100,'-.', T, -results.IT.oo_.irfs.TBY_IT_EPS_UC_RoW(1:40)'.*100,  T, -results.ES.oo_.irfs.TBY_ES_EPS_UC_RoW(1:40)'.*100,':', T,0*(1:40),':k','LineWidth',3);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
set(h(5),'linewidth',1);
title('Trade balance to GDP','FontSize',20)
legend(h([1 2 3 4]),{'DE','FR','IT','ES'},'FontSize',20,'Orientation','horizontal','Position',[0.265 0.01 0.5 0.05])
set(gcf, 'Position', get(0, 'Screensize'));
savefig([wd '\figures\' 'Saving_RoW'])
print('-f6', '-depsc', [wd '\figures\' 'Saving_RoW.eps'])
close(figure(6));

%%
% ****************************************
% Bond premium EA (Euro depreciation)
% ****************************************
figure(7),
subplot(3,3,1)
plot(-results.DE.oo_.irfs.LYOBS_DE_EPS_BW_EA(1:40)'.*25,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LYOBS_FR_EPS_BW_EA(1:40)'.*25,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LYOBS_IT_EPS_BW_EA(1:40)'.*25, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LYOBS_ES_EPS_BW_EA(1:40)'.*25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
% legend('DE','FR','IT','ES')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP','FontSize',20)

subplot(3,3,2)
plot(-results.DE.oo_.irfs.LC_DE_EPS_BW_EA(1:40)'.*25,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LC_FR_EPS_BW_EA(1:40)'.*25,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LC_IT_EPS_BW_EA(1:40)'.*25, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LC_ES_EPS_BW_EA(1:40)'.*25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Consumption','FontSize',20)

subplot(3,3,3)
plot(-results.DE.oo_.irfs.LI_DE_EPS_BW_EA(1:40)'.*25,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LI_FR_EPS_BW_EA(1:40)'.*25,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LI_IT_EPS_BW_EA(1:40)'.*25, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LI_ES_EPS_BW_EA(1:40)'.*25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Investment','FontSize',20)

subplot(3,3,4)
plot(-results.DE.oo_.irfs.LN_DE_EPS_BW_EA(1:40)'.*25,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LN_FR_EPS_BW_EA(1:40)'.*25,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LN_IT_EPS_BW_EA(1:40)'.*25, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LN_ES_EPS_BW_EA(1:40)'.*25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Hours','FontSize',20)

subplot(3,3,5)
plot(-results.DE.oo_.irfs.LWR_DE_EPS_BW_EA(1:40)'.*25,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LWR_FR_EPS_BW_EA(1:40)'.*25,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LWR_IT_EPS_BW_EA(1:40)'.*25, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LWR_ES_EPS_BW_EA(1:40)'.*25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real wages','FontSize',20)

subplot(3,3,6)
plot(-results.DE.oo_.irfs.R_DE_EPS_BW_EA(1:40)'.*25,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.R_FR_EPS_BW_EA(1:40)'.*25,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.R_IT_EPS_BW_EA(1:40)'.*25, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.R_ES_EPS_BW_EA(1:40)'.*25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real interest rate','FontSize',20)

subplot(3,3,7)
plot(-results.DE.oo_.irfs.PHIYOBS_DE_EPS_BW_EA(1:40)'.*25,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.PHIYOBS_FR_EPS_BW_EA(1:40)'.*25,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.PHIYOBS_IT_EPS_BW_EA(1:40)'.*25, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.PHIYOBS_ES_EPS_BW_EA(1:40)'.*25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('GDP Inflation','FontSize',20)

subplot(3,3,8)
plot(-results.DE.oo_.irfs.LRER_DE_EPS_BW_EA(1:40)'.*25,'--','LineWidth',3), hold on, plot(-results.FR.oo_.irfs.LRER_FR_EPS_BW_EA(1:40)'.*25,'-.','LineWidth',3), hold on, plot(-results.IT.oo_.irfs.LRER_IT_EPS_BW_EA(1:40)'.*25, 'LineWidth', 2), hold on, plot(-results.ES.oo_.irfs.LRER_ES_EPS_BW_EA(1:40)'.*25,':','LineWidth',3), hold on, plot(0*(1:40),':k','LineWidth',1);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
title('Real exchange rate','FontSize',20)

subplot(3,3,9)
h=plot(T,-results.DE.oo_.irfs.TBY_DE_EPS_BW_EA(1:40)'.*25,'--', T, -results.FR.oo_.irfs.TBY_FR_EPS_BW_EA(1:40)'.*25,'-.', T, -results.IT.oo_.irfs.TBY_IT_EPS_BW_EA(1:40)'.*25,  T, -results.ES.oo_.irfs.TBY_ES_EPS_BW_EA(1:40)'.*25,':', T,0*(1:40),':k','LineWidth',3);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
set(h(5),'linewidth',1);
title('Trade balance to GDP','FontSize',20);
legend(h([1 2 3 4]),{'DE','FR','IT','ES'},'FontSize',20,'Orientation','horizontal','Position',[0.265 0.01 0.5 0.05]);
set(gcf, 'Position', get(0, 'Screensize'));
savefig([wd '\figures\' 'Bond_premium'])
print('-f7', '-depsc', [wd '\figures\' 'Bond_premium.eps'])
close(figure(7));

