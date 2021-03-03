% clear all; close all; clc;
cd(workdirectory);
wd=pwd;
Countries = {'DE','FR','IT','ES'};
for c = 1:length(Countries)
    path = sprintf('%s', Countries{c});
    file1 = 'gemc_annual_fit_1.fig';

    Keep{1} = sprintf('GYOBSA_%s', Countries{c});
    Keep{2} = sprintf('PHIYOBSA_%s', Countries{c});
    Keep{3} = sprintf('GCA_%s', Countries{c});
    Keep{4} = sprintf('GIA_%s', Countries{c});
    Keep{5} = sprintf('GNA_%s', Countries{c});
    Keep{6} = sprintf('TBA_%s', Countries{c});
   
    Vec = 1:length(Keep);
    
    cTitles{1} = 'GDP growth'; %sprintf('GYOBS_%s(t-lag)', Countries{i});
    cTitles{2} = 'GDP deflator'; %sprintf('GC_%s(t-lag)', Countries{i});
    cTitles{3} = 'Consumption growth'; %sprintf('GI_%s(t-lag)', Countries{i});
    cTitles{4} = 'Investment growth'; %sprintf('PHIYOBS_%s(t-lag)', Countries{i});
    cTitles{5} = 'Growth of hours worked'; %sprintf('GN_%s(t-lag)', Countries{i});
    cTitles{6} = 'Trade balance'; %sprintf('GTBY_%s(t-lag)', Countries{i});
    
    
    % collecting the graphs
    fig1 = openfig(sprintf('%s\\%s', path, file1));
    for i=1:length(fig1.Children)
        ax = fig1.Children(i);
        if(ismember(ax.Title.String, Keep))
            ax_list(Vec(ismember(Keep, ax.Title.String))) = fig1.Children(i);
        end
    end
 
    % plot all together
    newfig = figure('Name', Countries{c}, 'units','normalized','outerposition',[0 0 1 1]);
    for i=1:6
       set(ax_list(i), 'parent', newfig);
       subplot(3,2,i,ax_list(i));
       title(cTitles{i}, 'FontSize', 12, 'Interpreter', 'latex');
    end
    close(fig1); 
    
    % Modify the graphs
    for i=1:6
        dataObjs = get(ax_list(i), 'Children');
        i
        % Find the lines to remove and treat
        for j=1:length(dataObjs)
            if( dataObjs(j).Color == [0 1 0] | dataObjs(j).Color == [0 1 1])% Green and light blue
                set(dataObjs(j), 'Visible', 'off'); % Remove
            end      
            if( dataObjs(j).Color == [0 0.4470 0.7410])% Blue
                set(dataObjs(j), 'LineWidth', 2);
                set(dataObjs(j), 'Color','black');
            end     
%             if( dataObjs(j).Color == [1 0 0])% red
%                 set(dataObjs(j), 'LineWidth', 2);
%                 set(dataObjs(j), 'Color', [0.5 0.5 0.5]);
%             end            
        end
    end
    %saveas(newfig,sprintf('Note_Annual_Fit_%s.png', Countries{c}));
%     savefig([wd '\figures\'  'Annual_Fit_%s.fig', Countries{c}]);
    cd('figures')
    saveas(newfig,sprintf('Annual_Fit_%s.fig', Countries{c}));
    saveas(newfig,sprintf('Annual_Fit_%s.eps', Countries{c}),'epsc');
    cd(wd)
    %print(fig_h, '-depsc','-r700',file_saveas)
    close(newfig);
end