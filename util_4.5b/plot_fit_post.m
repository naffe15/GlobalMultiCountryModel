function [rmse, lgobs_, r2, vv] = plot_fit_post(varargin);
%function [rmse, lgobs_] = plot_fit(varargin);
global options_ bayestopt_ M_ oo_

dirname = CheckPath('graphs',M_.dname);


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
x = T(fobs+istart-1:fobs+options_.nobs-1);
for ifig = 1:ceil(size(lgobs_,1)/9),
  h1 = dyn_figure(options_.nodisplay,'Name',['Bayesian DSGE 1-step ahead prediction']);
  h2 = dyn_figure(options_.nodisplay,'Name',['Bayesian DSGE Innovations']);
  for j=1+9*(ifig-1):min(9*ifig,size(lgobs_,1)),
    set(0,'CurrentFigure',h1)
    subplot(3,3,j-9*(ifig-1)),
    set(gca,'box','on')
    vj = deblank(lgobs_(j,:));
    i = strmatch(vj,lgy_,'exact');
    hp=patch([x; x(end:-1:1)], [oo_.FilteredVariables.HPDinf.(vj)(istart-1:end-1); oo_.FilteredVariables.HPDsup.(vj)(end-1:-1:istart-1)],[0.85 0.85 0.85],'EdgeColor',[0.85 0.85 0.85]);
    hold on,
    plot(x, oo_.FilteredVariables.Mean.(vj)(istart-1:end-1),'b','linewidth',1.5)
    plot(x, oo_.UpdatedVariables.Mean.(vj)(istart:end),'.-k')
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

        if options_.TeX,
            title(vj,'interpreter','none')
        else
            title(deblank(texname(i,:)),'interpreter','tex')
        end

%     hh=get(gca,'children');
%     set(hh(end-1),'linestyle','-','color',[0.7 0.7 0.7])
%     set(hh(end),'linestyle','-','color',[0 0 0])
%     hold on, plot([T(fobs)-1 T(fobs+options_.nobs-1)+1],[get_mean(vj) get_mean(vj)],'k:')
    plot(x, oo_.Constant.Mean.(vj)(istart:end),':k')
    set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
        dylim=max((yM-ym)*0.05,1.e-10);
        set(gca,'ylim',[ym-dylim yM+dylim])
    hold off,
    %disp([lgobs_(j,:), sprintf('%15.5g',[rmse(j)'])])
    figure(h2)
    subplot(3,3,j-9*(ifig-1)),
    vv(:,j) = oo_.UpdatedVariables.Mean.(vj)(istart:end)-oo_.FilteredVariables.Mean.(vj)(istart-1:options_.nobs-1);
    rmse(j,1) = sqrt(mean(vv(:,j).^2));
    r2(j,1) = 1-sum(vv(:,j).^2)/sum((oo_.UpdatedVariables.Mean.(vj)(istart:end)-oo_.Constant.Mean.(vj)(istart:end)).^2);    
%     r2(j,1) = 1-sum(vv(:,j).^2)/sum((oo_.UpdatedVariables.Mean.(vj)(istart:end)-trend(j,istart:end)').^2);    
    plot(x, vv(:,j),'k')
        if options_.TeX,
            title(vj,'interpreter','none')
        else
            title(deblank(texname(i,:)),'interpreter','tex')
        end

    set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
  end
  figure(h1)
  dyn_saveas(h1,[dirname,filesep,fname_,'_Bayesian_Fit',int2str(ifig)],options_.nodisplay,options_.graph_format)
  dyn_saveas(h2,[dirname,filesep,fname_,'_Bayesian_Innovations',int2str(ifig)],options_.nodisplay,options_.graph_format)
end


disp('Posterior FIT')
disp('             RMSE           R2')
% vname=str2mat(varargin{:});
for j=1:size(lgobs_,1),
    %     iv = strmatch(deblank(lgobs_(j,:)), lgobs_,'exact');
    disp([lgobs_(j,:), sprintf('%15.5g',[rmse(j)' r2(j)])])
end,
