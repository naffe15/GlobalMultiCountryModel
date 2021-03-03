function [rmse, lgobs_, r2, vv] = plot_fit(varargin);
%function [rmse, lgobs_] = plot_fit(varargin);
global options_ bayestopt_ M_ oo_

ys_ = oo_.steady_state;
lgy_ = M_.endo_names;
fname_ = M_.fname;
texname = M_.endo_names_tex;

if exist(options_.datafile)
    instr = options_.datafile;
else
    instr = ['load ' options_.datafile];
end
try
    eval(instr);
    if ~exist('T','var'),
        temp = eval(deblank(options_.varobs{1}));
        T=[1:length(temp)]';
        clear temp;
    end
catch
    T=[1:options_.nobs];
end
istart = max(2,options_.presample+1);

if options_.TeX,
    fidTeX = fopen([M_.fname '_Fit_Plots.TeX'],'w+');
end
if isempty(varargin),
    lgobs_= options_.varobs;
    mfys=bayestopt_.mfys;
else
    mfys=[];
    for j=1:length(varargin),
        dum = strmatch(varargin{j},lgy_,'exact');
        mfys = [mfys dum];
        if j==1,
            lgobs_ = varargin{j};
        else
            lgobs_ = str2mat(lgobs_,varargin{j});
        end
    end
end

nobs = size(lgobs_,1);
fobs = options_.first_obs;
if options_.loglinear == 1
    constant = log(ys_(mfys));
else
    constant = ys_(mfys);
end

trend = constant*ones(1,options_.nobs);

%disp(' ')
%disp(' ')
%disp('            RMSE')
numberofvars=nobs;

for ifig = 1:ceil(size(lgobs_,1)/9),
    
    
%     if options_.nodisplay,
%         h1 = figure('Name',['DSGE 1-step ahead prediction'],'visible','off');
%         h2 = figure('Name',['DSGE Innovations'],'visible','off');
%     else
        h1 = dyn_figure(options_.nodisplay,'Name',['DSGE 1-step ahead prediction']);
        h2 = dyn_figure(options_.nodisplay,'Name',['DSGE Innovations']);
%     end
    for j=1+9*(ifig-1):min(9*ifig,size(lgobs_,1)),
%         if options_.nodisplay
%             figure(h1,'visible','off')
%         else
        set(0,'CurrentFigure',h1)
%         end
        subplot(3,3,j-9*(ifig-1)),
        vj = deblank(lgobs_(j,:));
        i = strmatch(vj,lgy_,'exact');
        posi(j)=i;
        if ~isempty(texname),
            if strmatch('\log',texname(i,:));
                texname(i,1:end-1)=texname(i,2:end);
            end
        end
        offset=0;
        if exist(varargin{j},'var')
            eval(['plot(T(fobs+istart-1:fobs+options_.nobs-1),',varargin{j},'(fobs+istart-1:fobs+options_.nobs-1),''k.'')'])
            set(gco,'color',[0.7 0.7 0.7])
            hold on,
            offset=1;
        end
        eval(['plot(T(fobs+istart-1:fobs+options_.nobs-1), [oo_.FilteredVariables.', ...
            vj,'(istart-1:options_.nobs-1) oo_.UpdatedVariables.',...
            vj,'(istart:options_.nobs)])'])
        hh=get(gca,'children');
        ym = inf;
        yM = -inf;
        for ih=1:length(hh),
            try
            ym = min(ym,min(get(hh(ih),'ydata')));
            yM = max(yM,max(get(hh(ih),'ydata')));
            catch
            end
        end
        set(hh(end-offset-1),'linestyle','-','color',[0.7 0.7 0.7])
        set(hh(end-offset),'linestyle','-','color',[0 0 0]),
        if offset,
            set(hh(end),'color',[0.7 0.7 0.7]),
        end
        hold on, plot([T(fobs)-1 T(fobs+options_.nobs-1)+1],[get_mean(vj) get_mean(vj)],'r--')
        set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
        dylim=max((yM-ym)*0.05,1.e-10);
        set(gca,'ylim',[ym-dylim yM+dylim])
        if options_.TeX,
            title(vj,'interpreter','none')
        else
            title(deblank(texname(i,:)),'interpreter','tex')
        end
        eval(['rmse(j) = sqrt(mean((oo_.UpdatedVariables.',vj,'(istart:options_.nobs)-oo_.FilteredVariables.',vj,'(istart-1:options_.nobs-1)).^2));'])
        eval(['vv(:,j) = (oo_.UpdatedVariables.',vj,'(istart:options_.nobs)-oo_.FilteredVariables.',vj,'(istart-1:options_.nobs-1));'])
        eval(['r2(:,j) = 1-sum(vv(:,j).^2)/sum((oo_.UpdatedVariables.',vj,'(istart:options_.nobs)-trend(j,istart:options_.nobs)'').^2);'])
%         eval(['r2(:,j) = 1-cov(vv(:,j))/cov(oo_.UpdatedVariables.',vj,'(istart:end));'])
        %disp([lgobs_(j,:), sprintf('%15.5g',[rmse(j)'])])
%         if options_.nodisplay
%             figure(h2,'visible','off')
%         else
        set(0,'CurrentFigure',h2)
%         end
        subplot(3,3,j-9*(ifig-1)),
        plot(T(fobs+istart-1:fobs+options_.nobs-1), vv(:,j),'k')
        if options_.TeX,
            title(vj,'interpreter','none')
        else
            title(deblank(texname(i,:)),'interpreter','tex')
        end
        hold on, plot([T(fobs)-1 T(fobs+options_.nobs-1)+1],[0 0],'r--')
        set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
    end
    set(0,'CurrentFigure',h1)
%     end
%     eval(['print -depsc2 ' fname_,'_Fit',int2str(ifig)]);
%     eval(['print -dpdf ' fname_,'_Fit',int2str(ifig)]);
%     if options_.nodisplay
%         set(h1,'visible','on')
%     end
    dyn_saveas(h1,[fname_,'_Fit',int2str(ifig)],options_.nodisplay,options_.graph_format)
%     if options_.nodisplay
%         close(h1)
%     end
    set(0,'CurrentFigure',h2)
%     eval(['print -depsc2 ' fname_,'_Innovations',int2str(ifig)]);
%     eval(['print -dpdf ' fname_,'_Innovations',int2str(ifig)]);  
%     if options_.nodisplay
%         set(h2,'visible','on')
%     end
    dyn_saveas(h2,[fname_,'_Innovations',int2str(ifig)],options_.nodisplay,options_.graph_format)
%     if options_.nodisplay
%         close(h2)
%     end
    
    
    if options_.TeX,
        l1=9*(ifig-1)+1;
        texname0=texname(posi,:);
        %         fidTeX = fopen([M_.fname '_Fit_Plots.TeX'],'w+');
        %         files=dir('*_Fit*.eps');
%         fidTeX = fopen([M_.fname '_Fit_Plots.TeX'],'a+');
        tempfilename=[fname_,'_Fit',int2str(ifig)];
        [a,b,c]=fileparts(tempfilename);
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        
        for kk=1:min(numberofvars,9)
            
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(varargin{l1}),['$' deblank(texname0(l1,:)) '$']);
            l1=l1+1;
        end
        numberofvars=numberofvars-9;
        
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,['\\includegraphics[width=0.80\\textwidth] {' b '} \n']);
        fprintf(fidTeX,'\\caption{Fit Plot}');
        fprintf(fidTeX,'\\label{Fig: Fit Plot:%s}\n',int2str(ifig));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
    end
    
end

disp('             RMSE           R2')
% vname=str2mat(varargin{:});
for j=1:size(lgobs_,1),
    %     iv = strmatch(deblank(lgobs_(j,:)), lgobs_,'exact');
    disp([lgobs_(j,:), sprintf('%15.5g',[rmse(j)' r2(j)])])
end,

if options_.TeX,
    fclose(fidTeX);
    print_rmses_latex_table(varargin,rmse,r2);
end
    

