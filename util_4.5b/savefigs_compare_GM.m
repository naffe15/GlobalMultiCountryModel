
function savefigs_compare_GM(folders, countries1, vrbls_frcst, var_irf, shocks, vars_shock_decomp, shock_group)
mkdir compare_GM
wd=[pwd filesep 'compare_GM'];
% cd(wd)
mkdir(['compare_GM' filesep 'irfs'])
mkdir(['compare_GM' filesep 'posteriors'])
mkdir(['compare_GM' filesep 'frcst'])
% mkdir compare_GM' filesep 'zeps'])
mkdir(['compare_GM' filesep 'sDecomp'])
mkdir(['compare_GM' filesep 'fit'])

fnames = {
    'new',
    'base'};
countries = cellstr(countries1);
%
%% IRFS

for c = 1:length(countries)
    Cs    = cell(1, length(folders));
    Cs(:) = {countries{c}};
    for i = 1: length(shocks)
        figure, compare_GM(folders, var_irf,'irfs',shocks{i},Cs) ;
        h =  findobj('type','figure');
        for f = 1: length(h)
            figure(f);
            set(gcf, 'name', strcat(shocks{i},num2str(f),countries{c}));
            saveas(figure(f),[wd strcat('/irfs/',shocks{i},num2str(f),countries{c},'.fig')]);
        end
        close all
    end
end

%% POSTERIORS

if (exist([folders{1}  filesep 'gemc' filesep 'Output' filesep 'gemc_PriorsAndPosteriors1.fig' ]) == 2) && (exist([folders{2}  filesep 'gemc' filesep 'Output' filesep 'gemc_PriorsAndPosteriors1.fig' ])==2)
    
    figure, compare_GM(folders,'all','post_dens');
    h =  findobj('type','figure');
    for f = 1: length(h)
        figure(f);
        set(gcf, 'name', strcat('posterior',num2str(f)));
        saveas(figure(f),[wd strcat('/posteriors/','posterior',num2str(f),'.fig')]);
    end
    close all
end
%% Forecast
cd(folders{1})
vrbls_frcst =strcat(vrbls_frcst,countries);
figure, compare_GM(folders,vrbls_frcst,'frcst', 8);
cd(wd)
h =  findobj('type','figure');
for f = 1: length(h)
    figure(f);
    set(gcf, 'name', strcat('frcst',num2str(f)));
    saveas(figure(f),[wd strcat('/frcst/','frcst',num2str(f),'.fig')]);
end
close all
%
% cd(folders{1})
% figure, compare_GM(folders,vrbls_frcst1,'frcst', 8);
% cd(wd)
% h =  findobj('type','figure');
% for f = 1: length(h)
%     figure(f);
%     set(gcf, 'name', strcat('frcst1',num2str(f)));
%     saveas(figure(f),[pwd strcat('/frcst/','frcst1',num2str(f),'.fig')]);
% end

close all
% %% ZEPS
% % loading the results
% for f = 1:length(folders)
%     s{f}=strcat('f',num2str(f));
%     folders1{f}= strcat(folders{f},'\gemc_results.mat');
%     results.(s{f})=load(folders1{f});
%     lgnd{f} = erase(folders{f},'Z:\Global_Estimated_Model\GM\');
% end
% % plotting figures
% for i = 1:length(zeps)
%     for f =1:length(folders)
%         figure(i), plot([results.(s{f}).oo_.SmoothedVariables.(zeps{i})]);
%         hold on
%     end
%     legend(lgnd)
%     figure(i);
%     set(gcf, 'name', strcat(zeps{i}));
%     saveas(figure(i),[pwd strcat('/zeps/',zeps{i},'.fig')]);
% end
% close all
%%

% copy shock decomp
if (exist([folders{1}  filesep 'gemc' filesep 'graphs' filesep 'gemc_shock_decomposition_' 'GYOBS_' countries{c} '_yoy_group_' shock_group '.fig']) == 2) && (exist([folders{2}  filesep 'gemc' filesep 'graphs' filesep 'gemc_shock_decomposition_' 'GYOBS_' countries{c} '_yoy_group_' shock_group '.fig']) == 2)
    sddir = strcat(wd, [filesep 'sDecomp']);
    for c = 1:length(countries)
        for f = 1:length(folders)
            for i= 1:length(vars_shock_decomp)
                cd(strcat(folders{f}, [filesep 'gemc' filesep 'graphs']))
                if  strmatch(vars_shock_decomp{i},'TBY')==1
                    copyfile(strcat('gemc_shock_decomposition_',vars_shock_decomp{i},'_',countries{c},'_qoq_group_',shock_group,'.fig'),sddir);
                    cd(sddir)
                    copyfile(strcat('gemc_shock_decomposition_',vars_shock_decomp{i},'_',countries{c},'_qoq_group_',shock_group,'.fig'), strcat(vars_shock_decomp{i},'_',countries{c},'_',fnames{f},'.fig'));
                    delete(strcat('gemc_shock_decomposition_',vars_shock_decomp{i},'_',countries{c},'_qoq_group_',shock_group,'.fig'))
                else
                    copyfile(strcat('gemc_shock_decomposition_',vars_shock_decomp{i},'_',countries{c},'_yoy_group_',shock_group,'.fig'),sddir);
                    cd(sddir)
                    copyfile(strcat('gemc_shock_decomposition_',vars_shock_decomp{i},'_',countries{c},'_yoy_group_',shock_group,'.fig'), strcat(vars_shock_decomp{i},'_',countries{c},'_',fnames{f},'.fig'));
                    delete(strcat('gemc_shock_decomposition_',vars_shock_decomp{i},'_',countries{c},'_yoy_group_',shock_group,'.fig'))
                    
                end
                
                openfig(strcat(vars_shock_decomp{i},'_',countries{c},'_',fnames{f},'.fig'));
                set(gcf, 'name', strcat(vars_shock_decomp{i},'_',countries{c},'_',fnames{f},'.fig'));
                h=findobj('type','figure');
                saveas(figure(length(h)),[sddir strcat( filesep ,vars_shock_decomp{i},'_',countries{c},'_',fnames{f},'.fig')]);
                close all
            end
        end
    end
    cd(wd)
else
end
%%
% copy annual fit
sddir = strcat(wd, [filesep 'fit']);
for f = 1:length(folders)
    cd(folders{f})
    for i=1:2
        copyfile(fullfile(folders{f}, strcat('gemc_annual_fit','_',num2str(i),'.fig')),sddir);
        cd(sddir)
        openfig(strcat('gemc_annual_fit','_',num2str(i),'.fig'));
        %     h =  findobj('type','figure');
        %     figure(f);
        set(gcf, 'name', strcat('fit_',num2str(i),'_',fnames{f}));
        saveas(figure(i),[sddir strcat( filesep, 'fit_',num2str(i),'_',fnames{f},'.fig')]);
        delete(strcat('gemc_annual_fit','_',num2str(i),'.fig'));
        
    end
    close all
end

cd(folders{1})
%% copy r2a
% loading the results
for f = 1:length(folders)
    s{f}=strcat('f',num2str(f));
    folders1{f}= strcat(folders{f}, [filesep 'gemc_results.mat']);
    results.(s{f})=load(folders1{f});
    folders2{f}= strcat(folders{f}, [filesep 'gemc_results_cfrcst.mat']);
    results2.(s{f})=load(folders2{f},'annual_vname','r2a');
end

%saving the r2a

for f =1:length(folders)
    box{1,2*f}=fnames{f};
    if isfield(results2.(s{f}),'r2a')
        box(2:length(results2.(s{f}).annual_vname)+1,2*f)=results2.(s{f}).annual_vname;
        box(2:length(results2.(s{f}).annual_vname)+1,2*f+1)=num2cell(results2.(s{f}).r2a(:,1));
    else
        box(2:length(results.(s{f}).oo_.jrc.fit_annual.vname)+1,2*f)=results.(s{f}).oo_.jrc.fit_annual.vname;
        box(2:length(results.(s{f}).oo_.jrc.fit_annual.vname)+1,2*f+1)=num2cell(results.(s{f}).oo_.jrc.fit_annual.r2a(:,1));
    end
end
if ismac
xlswrite_MACOS([wd filesep 'r2a.xls'], box)
else
xlswrite([wd filesep 'r2a.xls'], box)
end
%% copy forecastbox
if (exist([folders{1} 'ForecastBox.xls' ]) == 2) && (exist([folders{2} 'ForecastBox.xls' ])==2)
    
    for f = 1:length(folders)
%         cd(folders{f});
        [A,B] = xlsfinfo([folders{f}  filesep 'ForecastBox.xls' ]);
        if any(strcmp(B, ['GYOBSA_' countries{1} '_Cond(No assumption)']))
            [numTr,txtTr,rawTr]=xlsread('ForecastBox.xls',['GYOBSA_' countries{1} '_Cond(No assumption)'],'A1:D20');
        else
            [numTr,txtTr,rawTr]=xlsread('ForecastBox.xls','Conditional(No assumption)','A1:D20');
            
        end
        s{f}=strcat('f',num2str(f));
        frcstbox.(s{f}).numTr = numTr;
        frcstbox.(s{f}).txtTr = txtTr;
    end
    cd(wd)
    for f =1:length(folders)
        fbox{1,5*(f-1)+2}=fnames{f};
        fbox(2:3,5*(f-1)+2:5*(f-1)+1+size(frcstbox.(s{f}).txtTr,2))=frcstbox.(s{f}).txtTr(1:2,:);
        
        fbox(2:length(frcstbox.(s{f}).txtTr)+1,5*(f-1)+2)=frcstbox.(s{f}).txtTr(:,1);
        fbox(5:length(frcstbox.(s{f}).numTr)+4,5*(f-1)+3:5*(f-1)+2+size(frcstbox.(s{f}).numTr,2))=num2cell(frcstbox.(s{f}).numTr);
    end
if ismac
    xlswrite_MACOS([wd filesep 'FBox.xls'], fbox)
else
    xlswrite([wd filesep 'FBox.xls'], fbox)
end
%     cd(folders{1})
end
end
