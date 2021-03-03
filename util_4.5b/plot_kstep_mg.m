function [rmseK, lgobs_] = plot_kstep(k_, varargin);
global options_ bayestopt_ M_ oo_

ys_ = oo_.steady_state;
lgy_ = M_.endo_names;
fname_ = M_.fname;

if isempty(k_), k_=options_.nk; end
if ischar(k_), 
    k_=eval(k_); 
end

if exist(options_.datafile)
    instr = options_.datafile;
else
    instr = ['load ' options_.datafile];
end
eval(instr);

if ~exist('T','var'), 
  temp = eval(deblank(options_.varobs{1}));
  T=[1:length(temp)]'; 
  clear temp;
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

istart=k_+options_.presample;
for ifig = 1:ceil(size(lgobs_,1)/9),
    dyn_figure(options_.nodisplay,'Name',['DSGE ',int2str(k_),'-step ahead prediction']);
    for j=1+9*(ifig-1):min(9*ifig,size(lgobs_,1)),
        subplot(3,3,j-9*(ifig-1)),
        vj=deblank(lgobs_(j,:));
%         i = strmatch(vj,lgy_(oo_.dr.order_var,:),'exact');
        i = strmatch(vj,lgy_,'exact');
        kk_ = find(options_.filter_step_ahead==k_);
        eval(['plot(T(fobs+istart-1:fobs+options_.nobs-1), [squeeze(oo_.FilteredVariablesKStepAhead', ...
                '(kk_,i,istart:end-options_.nk)) oo_.UpdatedVariables.',...
                vj,'(istart:end)])'])
%         eval(['plot(T(fobs+istart-1:fobs+options_.nobs-1), [oo_.FilteredVariables.', ...
%                 vj,'(istart:options_.nobs)+trend(j,istart:end)'' oo_.UpdatedVariables.',...
%                 vj,'(istart:end)+trend(j,istart:end)''])'])
        title(vj,'interpreter','none')
        hh=get(gca,'children');
        set(hh(1),'linestyle',':','color',[0 0 0])
        set(hh(2),'linestyle','-','color',[0 0 0])
        set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
        eval(['rmseK(j) = sqrt(mean((oo_.UpdatedVariables.',vj,'(istart:options_.nobs)-squeeze(oo_.FilteredVariablesKStepAhead(kk_,i,istart:end-options_.nk))).^2));'])
        %disp([lgobs_(j,:), sprintf('%15.5g',[rmseK(j)'])])
    end
    dyn_saveas(gcf,[fname_,'_',int2str(k_),'_pred',int2str(ifig)],options_.nodisplay,options_.graph_format)
end

if nargout>0,
disp(' ')
disp(' ')
disp('            k-step RMSE')
for j=1:size(lgobs_,1),
        disp([lgobs_(j,:), sprintf('%15.5g',[rmseK(j)'])])
end
end
