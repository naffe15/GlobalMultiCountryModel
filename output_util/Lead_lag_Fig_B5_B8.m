% clear all; close all; clc;

cd(workdirectory);
wd=pwd;

currpath{1}='DE';
currpath{2}='FR';
currpath{3}='IT';
currpath{4}='ES';

Countries = {'DE', 'FR', 'IT', 'ES'};
range = 1:5;    %nlags = 2;
%range= 1:11;     %nlags = 5;

for i=1:4
%     figure
    % load data   
    cfile = sprintf('%s\\gemc_moments.mat',currpath{i});
    load(cfile);
    cfile
    
        
    ComponentsCCF{1} = sprintf('GYOBS_%s_GYOBS_%s', Countries{i}, Countries{i});
    ComponentsCCF{2} = sprintf('GYOBS_%s_GC_%s', Countries{i}, Countries{i});
    ComponentsCCF{3} = sprintf('GYOBS_%s_GI_%s', Countries{i}, Countries{i});
    ComponentsCCF{4} = sprintf('GYOBS_%s_PHIYOBS_%s', Countries{i}, Countries{i});
    ComponentsCCF{5} = sprintf('GYOBS_%s_GN_%s', Countries{i}, Countries{i});
    ComponentsCCF{6} = sprintf('GYOBS_%s_GTBY_%s', Countries{i}, Countries{i});

    cTitles{1} = 'GDP growth'; %sprintf('GYOBS_%s(t-lag)', Countries{i});
    cTitles{2} = 'Consumption growth'; %sprintf('GC_%s(t-lag)', Countries{i});
    cTitles{3} = 'Investment growth'; %sprintf('GI_%s(t-lag)', Countries{i});
    cTitles{4} = 'GDP deflator'; %sprintf('PHIYOBS_%s(t-lag)', Countries{i});
    cTitles{5} = 'Hours growth'; %sprintf('GN_%s(t-lag)', Countries{i});
    cTitles{6} = '$\Delta$ Trade balance to GDP'; %sprintf('GTBY_%s(t-lag)', Countries{i});

    % plot
    newfig = figure('Name', Countries{i}, 'units','normalized','outerposition',[0 0 1 1]);
    for j=1:length(variables_)
        cData = data_moments.ccf.(ComponentsCCF{j});
        cModel = moments.ccf.(ComponentsCCF{j});
        
        subplot(3,2,j);
        stem(cData(range,1)+0.05, cData(range,2), '-b','filled', 'LineWidth', 2);
        hold on;
        stem(cModel(range,1)-0.05, cModel(range,2), '--sk','filled', 'LineWidth', 2);
        hold on;
        plot(cData(range,1) , cData(range,3) , '--r');
        hold on;
        plot(cData(range,1) , cData(range,4) , '--r');

        title(cTitles{j}, 'FontSize', 12, 'Interpreter', 'latex');
        grid on;
        xticks(-2:2)
        xlabel('lag', 'FontSize', 10);
        c_yticks = yticks;
        yticks(c_yticks(1):0.2:c_yticks(end));    
        set(gca, 'FontSize', 10)
    end
    lgd = legend('data', 'model');
    lgd.FontSize = 12;
    cd('figures')
    saveas(newfig,sprintf('Lead_lag_GDP_%s.fig', Countries{i}));
    saveas(newfig,sprintf('Lead_lag_GDP_%s.eps', Countries{i}),'epsc');
    close(newfig);
    cd(wd)
    % clear
    clear cfile ComponentsCCF cTitles cData cModel;
end