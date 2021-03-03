function x = get_shock_stderr_by_name_locally(exoname,M)
% function x = get_shock_stderr_by_name_locally(exoname)
% returns the value of a shock identified by its name
%  
% INPUTS:
%   exoname:  shock name
%
% OUTPUTS
%   x:      shock value
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2006-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.



i = strmatch(exoname,M.exo_names,'exact');

if isempty(i)
    error(sprintf('Can''t find shock %s', exoname))
end

x = sqrt(M.Sigma_e(i,i));
