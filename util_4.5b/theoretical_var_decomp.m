function tvd_matrix = theoretical_var_decomp(var_list_,names_vars,ex_names_,opts)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global M_ oo_

if nargin<4 || isempty(opts),
   opts = struct();
end

if isfield(opts,'flip_table')
    flip_table = opts.flip_table;
else
    flip_table = 0;   
end
if isfield(opts,'names_shock')
    names_shock = opts.names_shock;
else
    names_shock = ex_names_;   
end    
    
  
% Theoretical var decomp
tvd         = oo_.gamma_y{h};

% Shock groups
ngroups0    = size(ex_names_,1);

% Initialize matrix
tvd_matrix  = ones(size(var_list_,1),size(ex_names_,1)+1);


for i=1:ngroups0,
    clear index,
    for ii=1:size(ex_names_{i},2),
        indbuf = strmatch(ex_names_{i}{ii},M_.exo_names,'exact');
        if ~isempty(indbuf),
            index(ii) = indbuf;
        elseif ~isempty(ex_names_{i}{ii}),
            error(['Shock name ',ex_names_{i}{ii}, ' not found.' ]);
        end
        if size(ex_names_{i},2) > 1
            share_ngroup = sum(tvd(:,index)');
        else
            share_ngroup = tvd(:,index);
        end
    end    

    tvd_matrix(:,i) = share_ngroup;
end 
% Other Shocks
tvd_matrix = tvd_matrix*100;
tvd_matrix(:,ngroups0+1) =tvd_matrix(:,ngroups0+1) -  sum(tvd_matrix(:,1:ngroups0),2);
if flip_table == 1
    tvd_matrix      = tvd_matrix';
end

% %% Latex Output
% if flip_table == 1
%     tvd_matrix      = tvd_matrix';
%     rowLabels       = names_shock;
%     columnLabels    = names_vars;
% else
%     columnLabels    = names_shock;
%     rowLabels       = names_vars;
% end


%matrix2latex(tvd_matrix, 'out.tex', 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
%matrix2latex(tvd_matrix, 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
end

