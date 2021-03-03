
function y0 = get_irf(exo,varargin)
global M_ oo_

ys_ = [oo_.steady_state]*0;
y0=zeros(length(eval(['oo_.irfs.' varargin{1} '_' exo]))+1,length(varargin));

if isempty(varargin)
    varlist = cellstr(M_.endo_names(1:M_.orig_endo_nbr,:));
end

[i_var,nvar] = varlist_indices(char(varargin{:}),M_.endo_names);


for j=1:nvar,
%     mfys = strmatch(varargin{j},lgy_,'exact');
    y0(:,j)=[0; eval(['oo_.irfs.' varargin{j} '_' exo ])']+ys_(i_var(j));
end

