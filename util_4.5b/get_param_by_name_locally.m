function x = get_param_by_name_locally(pname,M)
% function x = get_param_by_name_locally(pname)
% returns the value of a parameter identified by its name
%
% INPUTS:
%   pname:  parameter name
%   M    :  structure 
%
% OUTPUTS
%   x:      parameter value
%
% SPECIAL REQUIREMENTS
%   none


i = strmatch(pname,M.param_names,'exact');

if isempty(i)
    error(sprintf('Can''t find parameter %s', pname))
end

x = M.params(i);
