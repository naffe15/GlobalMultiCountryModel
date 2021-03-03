%% load and save CU from oo_
% clear all;
cd(workdirectory)
wd=pwd;
mkdir([wd '\figures'])

countries={'DE', 'FR', 'IT', 'ES'};

b = countries;
for i= 1:length(countries)
 results.(b{i})=load([countries{i} '\gemc_results.mat']);
 results.(b{i}).dataobs=load([countries{i} '\dataobs.mat'], ['CUOBS' '_' countries{i} ]);
end


%%
set(groot,'defaultLineLineWidth',2)
T           = [1999:0.25:2016.75]';  
figure(10),
for c=1:length(countries)

subplot(2,2,c)
plot(T,results.(b{c}).oo_.SmoothedVariables.(['CU' '_' countries{c}])(1:72)',T,results.(b{c}).dataobs.(['CUOBS' '_' countries{c}])(16:87),'-.')
xlim([2000 2017])
title(b{c},'FontSize',9)
end
legend('Model', 'Data');
legend('Orientation','horizontal');
set(legend,'Position',[0.265 0.01 0.5 0.05]);
savefig([wd '\figures\' 'CU']);
print('-f10', '-depsc', [wd '\figures\' 'CU.eps']);
close(figure(10));



