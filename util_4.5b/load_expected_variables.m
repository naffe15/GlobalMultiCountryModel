function load_expected_variables()
% function load_expected_variables()
%
% INPUTS:
%
% OUTPUTS
%
% SPECIAL REQUIREMENTS
%   none

global M_ 

y0 = cellstr(M_.endo_names);
y1=get_expected(y0{:});
for j=1:M_.endo_nbr
    assignin('caller',['E_' deblank(M_.endo_names(j,:))],y1(:,j));
end
