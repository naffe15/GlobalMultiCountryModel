function load_param_values()
% function load_param_values()
%
% INPUTS:
%
% OUTPUTS
%
% SPECIAL REQUIREMENTS
%   none

global M_

for j=1:M_.param_nbr
    assignin('caller',deblank(M_.param_names(j,:)),M_.params(j));
end
