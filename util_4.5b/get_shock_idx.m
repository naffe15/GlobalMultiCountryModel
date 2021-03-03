function idx_shock = get_shock_idx(shockname,M)

idx_shock               = strmatch(shockname,M.exo_names);

end

