  function [z, shock_names, my_initial_date] = compare_plot_shock_decomp(M_,oo_,options_,varlist)
% function plot_shock_decomposition(M_,oo_,options_,varlist)
% Plots the results of shock_decomposition
%
% INPUTS
%    M_:          [structure]  Definition of the model
%    oo_:         [structure]  Storage of results
%    options_:    [structure]  Options
%    varlist:     [char]       List of variables
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2016-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

%options0 = evalin('base','options_');
%options0 = options_;
options_.nodisplay = options_.plot_shock_decomp.nodisplay;
options_.graph_format = options_.plot_shock_decomp.graph_format;

% indices of endogenous variables
exist_varlist = 1;
if size(varlist,1) == 0
    exist_varlist = 0;
    if size( M_.endo_names,1) >= M_.orig_endo_nbr
        varlist = M_.endo_names(1:M_.orig_endo_nbr,:);
    else
        varlist = M_.endo_names;
    end
end

if isfield(options_.plot_shock_decomp,'init2shocks') % private trap for uimenu calls
    init2shocks=options_.plot_shock_decomp.init2shocks;
else
    init2shocks=[];
end

if isfield(options_.shock_decomp,'i_var') && (M_.endo_nbr>=M_.orig_endo_nbr)
    M_.endo_names = M_.endo_names(options_.shock_decomp.i_var,:);
    M_.endo_names_tex = M_.endo_names_tex(options_.shock_decomp.i_var,:);
    M_.endo_nbr = length( options_.shock_decomp.i_var );
end
try
    [i_var,nvar,index_uniques] = varlist_indices(varlist,M_.endo_names);
    if ~isempty(init2shocks)
        varlist_indices(char(init2shocks{:,1}),M_.endo_names);
    end
catch ME
    if isfield(options_.shock_decomp,'i_var')
        warning('shock decomp results for some input variable was not stored: I recompute all decompositions')
        M_ = evalin('base','M_');
        bayestopt_ = evalin('base','bayestopt_');
        estim_params_ = evalin('base','estim_params_');
        var_list_ = char();
        disp('recomputing shock decomposition ...')
        [oo_,M_]= shock_decomposition(M_,oo_,options_,var_list_,bayestopt_,estim_params_);
        if isfield(oo_,'realtime_shock_decomposition') || options_.plot_shock_decomp.realtime
            disp('recomputing realtime shock decomposition ...')
            oo_ = realtime_shock_decomposition(M_,oo_,options_,var_list_,bayestopt_,estim_params_);
        end
        if isfield(oo_,'initval_decomposition')
            disp('recomputing initval shock decomposition ...')
            oo_ = initial_condition_decomposition(M_,oo_,options_,0,bayestopt_,estim_params_);   
        end
        assignin('base','oo_',oo_)
        [i_var,nvar,index_uniques] = varlist_indices(varlist,M_.endo_names);
        options_.shock_decomp = rmfield(options_.shock_decomp,'i_var');
       % options0.shock_decomp = rmfield(options0.shock_decomp,'i_var');
    else
        rethrow(ME)
    end
end
varlist=varlist(index_uniques,:);

% if ~isfield(options0.shock_decomp,'i_var') && exist_varlist
%     if ~isfield(options0.plot_shock_decomp,'i_var')
%         options0.plot_shock_decomp.i_var = i_var;
%     else
%         options0.plot_shock_decomp.i_var = union(i_var, options0.plot_shock_decomp.i_var);
%     end
%     assignin('base','options_',options0)
% end

type=options_.plot_shock_decomp.type;
if isequal(type, 'aoa') && isfield(options_.plot_shock_decomp,'q2a') && isstruct(options_.plot_shock_decomp.q2a)
    q2avec=options_.plot_shock_decomp.q2a;
    if nvar>1
        for jv = 1:nvar
            my_varlist = deblank(varlist(jv,:));
            indv = strcmp(my_varlist,{q2avec.qname});
            options_.plot_shock_decomp.q2a =  q2avec(indv);
            plot_shock_decomposition(M_,oo_,options_,my_varlist);
        end
        return
    else
        indv = strcmp(varlist,{q2avec.qname});
        options_.plot_shock_decomp.q2a =  q2avec(indv);
    end
end

% number of variables
endo_nbr = M_.endo_nbr;

% number of shocks
nshocks = M_.exo_nbr;
% type = '';
fig_name='';
% detail_plot=0;
% realtime_=0; % 0 is standard; 1 is realtime (pool/vintage); 2 is conditional (pool/vintage); 3 is forecast (pool/vintage)
% vintage_=0; % 0 pool realtime/conditional; int: forecast/conditional shock decompositions
% forecast_=0;
% steadystate=0;
% write_xls=0;

if isfield(options_.plot_shock_decomp,'cum_diff') % private trap for uimenu calls
    cum_diff_decomp=options_.plot_shock_decomp.cum_diff;
else
    cum_diff_decomp=0;
end
if isfield(options_.plot_shock_decomp,'diff') % private trap for uimenu calls
    differentiate_decomp=options_.plot_shock_decomp.diff;
else
    differentiate_decomp=0;
end
if cum_diff_decomp
    differentiate_decomp=1;
end    
if isfield(options_.plot_shock_decomp,'flip') % private trap for uimenu calls
    flip_decomp=options_.plot_shock_decomp.flip;
else
    flip_decomp=0;
end
if isfield(options_.plot_shock_decomp,'expand') % private trap for uimenu calls
    expand=options_.plot_shock_decomp.expand;
else
    expand=0;
    options_.plot_shock_decomp.expand=0;
end

options_.plot_shock_decomp.initval=0;
if ~isempty(options_.plot_shock_decomp.fig_name)
    fig_name=[' ' options_.plot_shock_decomp.fig_name];
    if length(fig_name)>=8 && strcmp(fig_name(end-6:end),'initval')
        options_.plot_shock_decomp.initval=1;
    end
end

detail_plot=options_.plot_shock_decomp.detail_plot;
realtime_= options_.plot_shock_decomp.realtime;
vintage_ = options_.plot_shock_decomp.vintage;
forecast_ = options_.shock_decomp.forecast;
steadystate = options_.plot_shock_decomp.steadystate;
write_xls = options_.plot_shock_decomp.write_xls;

if vintage_
    forecast_ = min(forecast_,options_.nobs-vintage_);
end

initial_date = options_.initial_date;
if isempty(initial_date)
    if isempty(type)
        % we assume annual model
        initial_date = dates('1Y');
    else
        % we assume the sample starts in Q1
        initial_date = dates('1Q1');
    end
end

if isfield(options_.plot_shock_decomp,'q2a') % private trap for aoa calls
    q2a=options_.plot_shock_decomp.q2a;
    if isstruct(q2a) && isempty(fieldnames(q2a))
        q2a=0;
    end
else
    q2a=0;
end

switch realtime_
    
    case 0
        if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
            z = oo_.shock_decomposition;
        end
        fig_name1=fig_name;
        
    case 1 % realtime
        if vintage_
            if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
                z = oo_.realtime_shock_decomposition.(['time_' int2str(vintage_)]);
            end
            fig_name1=[fig_name ' realtime (vintage ' char(initial_date+vintage_-1) ')'];
        else
            if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
                z = oo_.realtime_shock_decomposition.pool;
            end
            fig_name1=[fig_name ' realtime (rolling)'];
        end
        
    case 2 % conditional
        if vintage_
            if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
                z = oo_.realtime_conditional_shock_decomposition.(['time_' int2str(vintage_)]);
            end
            initial_date = initial_date+vintage_-1;
            fig_name1=[fig_name ' ' int2str(forecast_) '-step ahead conditional forecast (given ' char(initial_date) ')'];
        else
            if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
                z = oo_.conditional_shock_decomposition.pool;
            end
            fig_name1=[fig_name ' 1-step ahead conditional forecast (rolling)'];
        end
        
    case 3 % forecast
        if vintage_
            if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
                z = oo_.realtime_forecast_shock_decomposition.(['time_' int2str(vintage_)]);
            end
            initial_date = initial_date+vintage_-1;
            fig_name1=[fig_name ' ' int2str(forecast_) '-step ahead forecast (given ' char(initial_date) ')'];
        else
            if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
                z = oo_.realtime_forecast_shock_decomposition.pool;
            end
            fig_name1=[fig_name ' 1-step ahead forecast (rolling)'];
        end
end


if ~isempty(init2shocks) && ~expand
    n=size(init2shocks,1);
    for i=1:n
        j=strmatch(init2shocks{i,1},M_.endo_names,'exact');
        if ~isempty(init2shocks{i,2})
            jj=strmatch(init2shocks{i,2},M_.exo_names,'exact');
            z(:,jj,:)= z(:,jj,:) + oo_.initval_decomposition (:,j,:);
        else
            z(:,end,:)= z(:,end,:) - oo_.initval_decomposition (:,j,:);
        end
        z(:,end-1,:)= z(:,end-1,:) - oo_.initval_decomposition (:,j,:);
        
    end    
end    


if isfield(oo_.dr,'ys')
    steady_state = oo_.dr.ys;
else
    steady_state = oo_.steady_state;
end

if isequal(type,'aoa') && isstruct(q2a)
    if (expand && isfield(options_.plot_shock_decomp,'zfull'))
        za = options_.plot_shock_decomp.zfull;
        genda = size(za,3);
        endo_nbra = size(za,1);
        [za, shock_names, M_] = make_the_groups(za,genda,endo_nbra,nshocks,M_, options_);
    end
    if realtime_
        % take all dates where realtime is saved
        qqq=options_.initial_date+options_.shock_decomp.save_realtime(:)-1;
        % take the first Q4 of saved realtime
        t0=min(options_.shock_decomp.save_realtime(qqq.time(:,2)==4));
        if isempty(t0)
            error('the realtime decompositions are not stored in Q4! Please check your dates and settings.')
        end
        if ~isfield(q2a,'var_type') % private trap for aoa calls
            q2a.var_type=1;
        end
        if ~isfield(q2a,'islog') % private trap for aoa calls
            q2a.islog=0;
        end
        if ~isfield(q2a,'GYTREND0') % private trap for aoa calls
            q2a.GYTREND0=0;
        end
        if ~isfield(q2a,'aux') % private trap for aoa calls
            q2a.aux=0;
        end
        if ~isfield(q2a,'cumfix') % private trap for aoa calls
            q2a.cumfix=1;
        end
        if ~isfield(q2a,'plot') % private trap for aoa calls
            q2a.plot=1; % growth rate
        end
        
        %     if isstruct(q2a.aux) && ischar(q2a.aux.y)
        %         opts=options_;
        %         opts.plot_shock_decomp.type='qoq';
        %         [y_aux, steady_state_aux] = plot_shock_decomposition(M_,oo_,opts,q2a.aux.y);
        %         q2a.aux.y=y_aux;
        %         q2a.aux.yss=steady_state_aux;
        %     end
        if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
            if options_.plot_shock_decomp.interactive && ~isempty(options_.plot_shock_decomp.use_shock_groups)
                mygroup = options_.plot_shock_decomp.use_shock_groups;
                options_.plot_shock_decomp.use_shock_groups='';
                zafull = ...
                    annualized_shock_decomposition(oo_,M_, options_, i_var, t0, options_.nobs, realtime_, vintage_, steady_state,q2a);
                options_.plot_shock_decomp.use_shock_groups = mygroup;
            end
            [za, endo_names, endo_names_tex, steady_state, i_var, oo_] = ...
                annualized_shock_decomposition(oo_,M_, options_, i_var, t0, options_.nobs, realtime_, vintage_, steady_state,q2a);
        end
        %     if realtime_<2
        %         initial_date = initial_date1;
        %     else
        %         initial_date = initial_date0;
        %     end
    end
end



if ~expand
    fig_name = fig_name1;
end
if options_.plot_shock_decomp.use_shock_groups
    fig_name=[fig_name ' group ' options_.plot_shock_decomp.use_shock_groups];
    if (expand && isfield(options_.plot_shock_decomp,'zfull'))
        z = options_.plot_shock_decomp.zfull;
        endo_names = options_.plot_shock_decomp.endo_names;
        endo_names_tex = options_.plot_shock_decomp.endo_names_tex;
        steady_state=steady_state(i_var);
        i_var = 1;
        if ~(isequal(type,'aoa') && isstruct(q2a))
            gend = size(z,3);
            endo_nbr = size(z,1);
            [z, shock_names, M_] = make_the_groups(z,gend,endo_nbr,nshocks,M_,options_);
            M_.endo_names = endo_names;
            M_.endo_names_tex = endo_names_tex;
        else
            % here we know we only have one variable to handle
            [~, yssa, ~, gyssa] = ...
                quarterly2annual(repmat(steady_state,16,1),steady_state,q2a.GYTREND0,q2a.type,q2a.islog,q2a.aux);
            if q2a.plot==1
                steady_state = gyssa;
            else
                steady_state = yssa;
            end
        end
    else
        gend = size(z,3);
        zfull = z;
        [z, shock_names, M_] = make_the_groups(z,gend,endo_nbr,nshocks,M_,options_);
    end
else
    shock_names = M_.exo_names;
end

func = @(x) colorspace('RGB->Lab',x);
MAP = distinguishable_colors(size(z,2)-1,'w',func);
%         MAP = [MAP; MAP(end,:)];
MAP(end,:) = [0.7 0.7 0.7];
%         MAP = [MAP; [0.7 0.7 0.7]; [0.3 0.3 0.3]];

if isempty(options_.plot_shock_decomp.colormap)
    options_.plot_shock_decomp.colormap = MAP;
end

if differentiate_decomp
    z(:,:,2:end) =  z(:,:,2:end)-z(:,:,1:end-1);
    z(:,:,1) = nan;
    if cum_diff_decomp
        z(:,:,1) = 0;
        z=cumsum(z,3);
    end
end

switch type

  case '' % default

  case 'qoq'

  case 'yoy'
    z=z(:,:,1:end-3)+z(:,:,2:end-2)+z(:,:,3:end-1)+z(:,:,4:end);
    if ~isempty(initial_date)
        initial_date = initial_date+3;
    else
        initial_date = dates('1Q4');
    end
    steady_state = 4*steady_state;

  case 'aoa'

    if isempty(initial_date)
        t0=1; % we assume the sample starts Q1 of 1st year
        initial_date = dates('1Y');
    else
        initial_date0 = dates([int2str(initial_date.time(1)) 'Y']);
        if initial_date.time(2)==1  % the first year is full
            t0=1;
            initial_date1=initial_date0;
        else
            t0=(4-initial_date.time(2)+2); % 1st period of the 1st full year in sample
            initial_date1=initial_date0+1;
        end
    end
    if realtime_ == 0
        t0=t0+4-1; % we start in Q4 of the first full year
        end
    if isempty(options_.plot_shock_decomp.plot_init_date) && realtime_ == 0
        options_.plot_shock_decomp.plot_init_date=initial_date+t0;
    end
    if isstruct(q2a)
        if realtime_ == 0
            if ~isfield(q2a,'var_type') % private trap for aoa calls
                q2a.var_type=1;
            end
            if ~isfield(q2a,'islog') % private trap for aoa calls
                q2a.islog=0;
            end
            if ~isfield(q2a,'GYTREND0') % private trap for aoa calls
                q2a.GYTREND0=0;
            end
            if ~isfield(q2a,'aux') % private trap for aoa calls
                q2a.aux=0;
            end
            if ~isfield(q2a,'cumfix') % private trap for aoa calls
                q2a.cumfix=1;
            end
            if ~isfield(q2a,'plot') % private trap for aoa calls
                q2a.plot=1; % growth rate
            end

            if ~(expand && isfield(options_.plot_shock_decomp,'zfull'))
            if isstruct(q2a.aux) && ischar(q2a.aux.y)
                opts=options_;
                opts.plot_shock_decomp.type='qoq';
                [y_aux, steady_state_aux] = plot_shock_decomposition(M_,oo_,opts,q2a.aux.y);
                q2a.aux.y=y_aux;
                q2a.aux.yss=steady_state_aux;
            end
            [za, endo_names, endo_names_tex, steady_state, i_var, oo_] = ...
                annualized_shock_decomposition(z,M_, options_, i_var, t0, options_.nobs, realtime_, vintage_, steady_state,q2a);
                if options_.plot_shock_decomp.interactive && ~isempty(options_.plot_shock_decomp.use_shock_groups)
                    mygroup = options_.plot_shock_decomp.use_shock_groups;
                    options_.plot_shock_decomp.use_shock_groups='';
                    zafull = ...
                        annualized_shock_decomposition(z,M_, options_, i_var, t0, options_.nobs, realtime_, vintage_, steady_state,q2a);
                    options_.plot_shock_decomp.use_shock_groups = mygroup;
                end
            end
        end
        z = za;
            if options_.plot_shock_decomp.interactive && ~isempty(options_.plot_shock_decomp.use_shock_groups)
                zfull = zafull;
            end
        M_.endo_names = endo_names;
        M_.endo_names_tex = endo_names_tex;
        %     endo_nbr = size(z,1);
        if realtime_<2 || vintage_ == 0
            initial_date = initial_date1;
        else
            initial_date = initial_date0;
        end
    else
	    % this is for quarterly-annualized variables already defined in model, so we can just take Q4
        t0=4-initial_date.time(2)+1;
        initial_date = initial_date0;
        z=z(:,:,t0:4:end);
    end

    if ~isempty(options_.plot_shock_decomp.plot_init_date)
        options_.plot_shock_decomp.plot_init_date = dates([int2str(options_.plot_shock_decomp.plot_init_date.time(1)) 'Y']);
    end
    if ~isempty(options_.plot_shock_decomp.plot_end_date)
        options_.plot_shock_decomp.plot_end_date = dates([int2str(options_.plot_shock_decomp.plot_end_date.time(1)) 'Y']);
    end


  otherwise

    error('plot_shock_decomposition:: Wrong type')

end

if flip_decomp
    z = -z;
    steady_state = - steady_state;
end
if steadystate
    options_.plot_shock_decomp.steady_state=steady_state;
end



% here we crop data if needed
my_initial_date = initial_date;
a = 1;
b = size(z,3);
if ~isempty(options_.plot_shock_decomp.plot_init_date)
    my_initial_date = max(initial_date,options_.plot_shock_decomp.plot_init_date);
    a = find((initial_date:initial_date+b-1)==options_.plot_shock_decomp.plot_init_date);
end
if ~isempty(options_.plot_shock_decomp.plot_end_date)
    b = find((initial_date:initial_date+b-1)==options_.plot_shock_decomp.plot_end_date);
end
z = z(:,:,a:b);
% end crop data
if nargout
    z=z(i_var,:,:);
    steady_state = steady_state(i_var);
    my_initial_date = initial_date;
    return
end
options_.plot_shock_decomp.fig_name=fig_name;
options_.plot_shock_decomp.orig_varlist = varlist;
if options_.plot_shock_decomp.interactive && ~isempty(options_.plot_shock_decomp.use_shock_groups)
    options_.plot_shock_decomp.zfull = zfull;
end
if detail_plot
    graph_decomp_detail(z,shock_names,M_.endo_names,i_var,my_initial_date,M_,options_)
else
    compare_graph_decomp(z,shock_names,M_.endo_names,i_var,my_initial_date,M_,options_);
end

if write_xls 
    WriteShockDecomp2Excel(z,shock_names,M_.endo_names,i_var,my_initial_date,M_,options_,options_.plot_shock_decomp);
end

end

function [z, shock_names, M_] = make_the_groups(z,gend,endo_nbr,nshocks,M_,options_)

shock_groups = M_.shock_groups.(options_.plot_shock_decomp.use_shock_groups);
shock_ind = fieldnames(shock_groups);
ngroups = length(shock_ind);
shock_names = shock_ind;
for i=1:ngroups
    shock_names{i} = (shock_groups.(shock_ind{i}).label);
end
zz = zeros(endo_nbr,ngroups+2,gend);
kcum=[];
for i=1:ngroups
    for j = shock_groups.(shock_ind{i}).shocks
        k = find(strcmp(j,cellstr(M_.exo_names)));
        zz(:,i,:) = zz(:,i,:) + z(:,k,:);
        z(:,k,:) = 0;
        kcum = [kcum k];
    end
end
zothers = sum(z(:,1:nshocks,:),2);
shock_groups.(['group' int2str(ngroups+1)]).label =  'Others';
shock_groups.(['group' int2str(ngroups+1)]).shocks =  cellstr(M_.exo_names(find(~ismember([1:M_.exo_nbr],kcum)),:))';
M_.shock_groups.(options_.plot_shock_decomp.use_shock_groups)=shock_groups;
if any(any(zothers))
    shock_names = [shock_names; {'Others + Initial Values'}];
end
zz(:,ngroups+1,:) = sum(z(:,1:nshocks+1,:),2);
zz(:,ngroups+2,:) = z(:,nshocks+2,:);
z = zz;

end
