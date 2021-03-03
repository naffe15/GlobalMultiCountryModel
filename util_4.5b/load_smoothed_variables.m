function load_smoothed_variables()
% function load_smoothed_variables()
%
% INPUTS:
%
% OUTPUTS
%
% SPECIAL REQUIREMENTS
%   none

global M_ oo_

for j=1:M_.endo_nbr
    assignin('caller',deblank(M_.endo_names(j,:)),oo_.SmoothedVariables.(deblank(M_.endo_names(j,:))));
end
