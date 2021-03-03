function xparam1 = calib2mode(bayestopt,M)

np = length(bayestopt.name);
xparam1=zeros(np,1);
for j=1:np,
    if strmatch(bayestopt.name{j},M.exo_names)
        xparam1(j) = get_shock_stderr_by_name(bayestopt.name{j});
    end
    if strmatch(bayestopt.name{j},M.param_names)
        xparam1(j) = get_param_by_name(bayestopt.name{j});
    end
end