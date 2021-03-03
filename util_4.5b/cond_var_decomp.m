function vd_tot = cond_var_decomp(var_list_,ex_names_,opts)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global M_ oo_ options_

if nargin<3 || isempty(opts),
   opts = struct();
end

if isfield(opts,'flip_table')
    flip_table = opts.flip_table;
else
    flip_table = 0;   
end

if isfield(opts,'names_shocks')
    names_shock = opts.names_shocks;
else
    names_shock = ex_names_;   
end

% Get vars
nvars    = size(var_list_,1);
% idx_vars = zeros(nvars,1);
% for j=1:nvars
%     idx_vars(j) = strmatch(var_list_(j,:),M_.endo_names);
% end
% 
if isfield(opts,'names_vars')
    names_vars = opts.names_vars;
else
    names_vars = M_.endo_names_tex(idx_vars,:);   
end

% Shock groups
ngroups0    = size(ex_names_,1);

% nvars       = size(cond_sd_mat,1);
% nhorizon    = size(cond_sd_mat,2);


%% Setup matrix
% cond_vdec_mat = oo_.conditional_variance_decomposition(idx_vars,:,:);
cond_vdec_mat = oo_.conditional_variance_decomposition;

% Initialize matrix
nhorizon =  size(options_.conditional_variance_decomposition,1);
tvd_matrix   =  ones(nvars,nhorizon,size(ex_names_,1)+1);
for i=1:ngroups0
    
    clear index,
    clear share_ngroup
    
    % Find shock
    % -----------
    for ii=1:size(ex_names_{i},2),
        indbuf = strmatch(ex_names_{i}{ii},M_.exo_names,'exact');
        if ~isempty(indbuf),
            index(ii) = indbuf;
        elseif ~isempty(ex_names_{i}{ii}),

            
        end
    end
        
    % Get shock group share
    % ----------------------
    if size(ex_names_{i},2) > 1
        share_ngroup = sum(cond_vdec_mat (:,:,index),3);
    else
        share_ngroup = cond_vdec_mat(:,:,index);
    end
    tvd_matrix(:,:,i) = share_ngroup;
end 
%Other Shocks
tvd_matrix_test =  1 - sum(tvd_matrix(:,:,1:ngroups0),3);
tvd_matrix(:,:,ngroups0+1) =tvd_matrix(:,:,ngroups0+1) - sum(tvd_matrix(:,:,1:ngroups0),3);
tvd_matrix = tvd_matrix*100;


%% Output        
vd_tot = reshape(tvd_matrix,size(tvd_matrix,1)*size(tvd_matrix,2),ngroups0+1);

% Print into LaTeX
print_cond_var_latex_table(vd_tot,ex_names_,names_shock,names_vars)
















% % tvd_tot = reshape(tvd_matrix,size(tvd_matrix,1)*size(tvd_matrix,2),10);
% % 
% % flip_table = 0;
% % if flip_table
% %     tvd_tot         = tvd_tot';
% %     rowLabels       = names_shock;
% %     columnLabels    = repmat(names_vars,1,nhorizon);  
% % else
% %     rowLabels       = repmat(names_vars,1,nhorizon);
% %     columnLabels    = names_shock;  
% % end
% % 
% % matrix = tvd_tot(:,:)
% % 
% % rowLabels = {'$g^y$','$g^{AY}$','$g^{XRD}$','$g^{AYMAX}$'} 
% % rowLabels = {'GDP growth',...
% %             'Productivity growth ',...
% %             'R\&D growth',...
% %             'Frontier growth'}
% % 
% % columnLabels = {'$g^y_1$','$g^{AY}_1$','$g^{XRD}_1$','$g^{AYMAX}_1$' '$g^y_4$','$g^{AY}_4$','$g^{XRD}_4$','$g^{AYMAX}_4$', '$g^y_\infty$','$g^{AY}_\infty$','$g^{XRD}_\infty$','$g^{AYMAX}_\infty$'};
% %   
% % matrix2latex(tvd_tot , 'out.tex', 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');  
% % matrix2latex(tvd_tot, 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

  
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % % Theoretical var decomp
% % tvd         = oo_.gamma_y{h};
% % 
% % % Shock groups
% % ngroups0    = size(ex_names_,1);
% % 
% % % Initialize matrix
% % tvd_matrix  = ones(size(var_list_,1),size(ex_names_,1)+1);
% % 
% % 
% % for i=1:ngroups0,
% %     clear index,
% %     for ii=1:size(ex_names_{i},2),
% %         indbuf = strmatch(ex_names_{i}{ii},M_.exo_names,'exact');
% %         if ~isempty(indbuf),
% %             index(ii) = indbuf;
% %         elseif ~isempty(ex_names_{i}{ii}),
% %             error(['Shock name ',ex_names_{i}{ii}, ' not found.' ]);
% %         end
% %         if size(ex_names_{i},2) > 1
% %             share_ngroup = sum(tvd(:,index)');
% %         else
% %             share_ngroup = tvd(:,index);
% %         end
% %     end    
% % 
% %     tvd_matrix(:,i) = share_ngroup;
% % end 
% % % Other Shocks
% % tvd_matrix = tvd_matrix*100;
% % tvd_matrix(:,ngroups0+1) =tvd_matrix(:,ngroups0+1) -  sum(tvd_matrix(:,1:ngroups0),2);
% % if flip_table == 1
% %     tvd_matrix      = tvd_matrix';
% % end
% % 
% % % %% Latex Output
% % % if flip_table == 1
% % %     tvd_matrix      = tvd_matrix';
% % %     rowLabels       = names_shock;
% % %     columnLabels    = names_vars;
% % % else
% % %     columnLabels    = names_shock;
% % %     rowLabels       = names_vars;
% % % end


%matrix2latex(tvd_matrix, 'out.tex', 'rowLabels', rowLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
%matrix2latex(tvd_matrix, 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
end

