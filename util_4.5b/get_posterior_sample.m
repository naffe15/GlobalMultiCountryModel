function [ oo, params, varlist ] = get_posterior_sample( type, varlist, M_, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   for type = 'param'
%   to get all params
%  [ oo, pp, varlist ] = get_posterior_sample(  'param', [], M_, bayestopt_ );
%   to get 1 or a list of params
%  [ oo, pp, varlist ] = get_posterior_sample(  'param', 'parname', M_, bayestopt_ );
%   to get posterior kernel values
%  [ oo, pp, varlist ] = get_posterior_sample(  'param', 'logpo2', M_, bayestopt_ );

    if isequal(type,'param')
        bayestopt_ = varargin{1};
    end
if size(varlist,1) == 0
    if isequal(type,'param')
        varlist = char(bayestopt_.name{:});
    else
        varlist = M_.endo_names(1:M_.orig_endo_nbr,:);
    end
end

    if ~isequal(type,'param')
        [i_var,nvar] = varlist_indices(varlist,M_.endo_names);
    end
switch type
    
    case 'irf'
        param_file = [M_.fname '_param_irf'];
        
    case 'param'
        
        param_file = [];
        DirectoryName = CheckPath('metropolis',M_.dname);
        a=dir([DirectoryName '/'  M_.fname '_mh_history*']);
        load([DirectoryName '/'  M_.fname '_mh_history_' int2str(length(a)-1)]);
        % load([M_.fname '/prior/definition.mat']);
        TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
        TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
        npar = size(varlist,1);
        
        for j=1:npar,
            jp = strmatch(deblank(varlist(j,:)),bayestopt_.name,'exact');
            params(:,j) = GetAllPosteriorDraws(jp,1,1,TotalNumberOfMhFiles,TotalNumberOfMhDraws);
        end
        
        oo=GetAllPosteriorDraws(0,1,1,TotalNumberOfMhFiles,TotalNumberOfMhDraws);
        return
        
    otherwise
        param_file = [M_.fname '_param'];
end

ifile = 1;
tmpfile = dir([M_.fname filesep 'metropolis' filesep param_file int2str(ifile) '.mat']);
if isempty(tmpfile)
    error(['no posterior sample for type: ' type])
    return
end

params=[];logpost=[];steady_state=[];
while ~ isempty(tmpfile)
    tmp_params = load([M_.fname filesep 'metropolis' filesep param_file int2str(ifile) '.mat']);
    params = [params; tmp_params.stock];
    if ~isequal(type, 'irf')
        logpost = [logpost; tmp_params.stock_logpo];
        steady_state = [steady_state; tmp_params.stock_ys];
    end
    ifile = ifile+1;
    tmpfile = dir([M_.fname filesep 'metropolis' filesep param_file int2str(ifile) '.mat']);
    
end


switch type
    
    case 'filter'
        
        type_file = [M_.fname '_filter_step_ahead'];
        
    case 'forecast'
        
        type_file = [M_.fname '_forc_mean'];
        
    case 'irf'
        
        type_file = [M_.fname '_IRF_DSGEs'];
        
    case 'smooth'
        
        type_file = [M_.fname '_smooth'];
        
    case 'update'
        
        type_file = [M_.fname '_update'];
        
        
end


jfile = 1;
tmpfile = dir([M_.fname filesep 'metropolis' filesep type_file int2str(jfile) '.mat']);
if isempty(tmpfile)
    error(['no posterior sample for type: ' type])
    return
end

switch type
    
    case 'irf'
        
        % Get index of shocks for requested IRFs
        if isempty(varargin)
            irf_shocks = M_.exo_names;
            irf_shocks_indx = M_.exo_names_orig_ord;
        else
            irf_shocks = varargin{1};
            irf_shocks_indx = zeros(1,size(irf_shocks,1));
            for i=1:size(irf_shocks,1)
                irf_shocks_indx(i) = find(strcmp(deblank(irf_shocks(i,:)), cellstr(M_.exo_names)));
            end
            irf_shocks_indx_unique=unique(irf_shocks_indx);
            irf_shocks_indx=irf_shocks_indx_unique;
        end
        oo=[];
        kdx = 0;
        while ~ isempty(tmpfile)
            load([M_.fname filesep 'metropolis' filesep type_file int2str(jfile) '.mat']);
            k = size(STOCK_IRF_DSGE,1);
            kk = (1:k)+kdx;
            ix=0;
            for i = irf_shocks_indx
                ix = ix+1;
                for j=1:length(i_var)
                    oo.(deblank(irf_shocks(ix,:)))(j,kk,:)=squeeze(STOCK_IRF_DSGE(:,i_var(j),i,:));
                end
            end
            kdx = kdx + size(STOCK_IRF_DSGE,1);
            
            jfile = jfile+1;
            tmpfile = dir([M_.fname filesep 'metropolis' filesep type_file int2str(jfile) '.mat']);
        end
        
    case 'filter'

        oo0=[];
        while ~ isempty(tmpfile)
            tmp_oo = load([M_.fname filesep 'metropolis' filesep type_file int2str(jfile) '.mat']);
            oo0 = cat(4, oo0, tmp_oo.stock(:,i_var,:,:));
            jfile = jfile+1;
            tmpfile = dir([M_.fname filesep 'metropolis' filesep type_file int2str(jfile) '.mat']);
            
        end
        for jv=1:nvar,
            oo.(deblank(varlist(jv,:)))= squeeze(oo0(:,jv,:,:));
        end
         
    otherwise
        
        oo=[];
        while ~ isempty(tmpfile)
            tmp_oo = load([M_.fname filesep 'metropolis' filesep type_file int2str(jfile) '.mat']);
            oo= cat(3, oo, tmp_oo.stock(i_var,:,:));
            jfile = jfile+1;
            tmpfile = dir([M_.fname filesep 'metropolis' filesep type_file int2str(jfile) '.mat']);
            
        end
        
        
        
end

end

