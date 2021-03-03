function [fcast] = plot_fcast(k_, varargin);
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
if ~exist('T'), 
  temp = eval(deblank(options_.varobs(1,:)));
  T=[1:length(temp)]'; 
  clear temp;
end

if isempty(varargin)
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

lst_obs = options_.nobs+options_.first_obs - 1;

%disp(' ')
%disp(' ')
%disp('            k-step RMSE')
for ifig = 1:ceil(size(lgobs_,1)/9),
    figure('Name',['DSGE ',int2str(k_),'-step forecasts']),
    for j=1+9*(ifig-1):min(9*ifig,size(lgobs_,1)),
        subplot(3,3,j-9*(ifig-1)),
        vj=deblank(lgobs_(j,:));
        i = strmatch(vj,lgy_,'exact');
        for l=1:k_
            eval(['vp(l)= oo_.KStepPredictions.', ...
                    vj,'(options_.nobs+l,l)+constant(j);'])
        end
        plot([T(lst_obs)+0.25:0.25:T(lst_obs)+k_/4], vp)
        hold on,
        if exist(vj,'var')
            
            eval(['vd= ',vj,';'])
            plot([T(lst_obs-8:end)], vd(lst_obs-8:end),'r')  
        else
%             eval(['vd= oo_.SmoothedVariables.',vj,'+constant(j);'])
            eval(['vd= oo_.SmoothedVariables.',vj,';'])
            plot(T(lst_obs-8:end), vd(end-8:end) ,'r')
            
        end
        title(vj,'interpreter','none')
        hh=get(gca,'children');
    end
    saveas(gcf,[fname_,'_',int2str(k_),'_fcast',int2str(ifig)])
    if options_.nodisplay,
        close(gcf)
    end
end
