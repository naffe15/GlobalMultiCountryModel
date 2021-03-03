function [endo_var_name, ind_aux_vars] = get_aux_var_parent_name(aux_var_name)
global M_ 

ind_aux = strmatch(aux_var_name,M_.endo_names);
ggg=[M_.aux_vars.endo_index];

ind_aux_vars = find(ggg==ind_aux);
if ~isempty(M_.aux_vars(ind_aux_vars).orig_index)
    endo_var_name = deblank(M_.endo_names(M_.aux_vars(ind_aux_vars).orig_index,:));
    if ~isempty(M_.aux_vars(ind_aux_vars).orig_lead_lag)
        endo_var_name = [endo_var_name '('  int2str(M_.aux_vars(ind_aux_vars).orig_lead_lag) ')'];
    end
else
    disp('Cannot deterine the parent of this aux var')
    disp('Check entry of M_')
    disp(['M_.aux_vars(' int2str(ind_aux_vars) ')'])
    disp('or the _dynamic.tex file')
    endo_var_name = [];
    return
end
