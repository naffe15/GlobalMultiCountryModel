function set_coefficients(xparam1)
  global estim_params_ M_
  
  nvx = estim_params_.nvx;
  ncx = estim_params_.ncx;
  nvn = estim_params_.nvn;
  ncn = estim_params_.ncn;
  np = estim_params_.np;
  
  if np
    offset = nvx+ncx+nvn+ncn;
    M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
  end
  
