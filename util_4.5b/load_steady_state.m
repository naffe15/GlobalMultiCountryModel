function load_steady_state()
% function load_steady_state()
%
% INPUTS:
%
% OUTPUTS
%
% SPECIAL REQUIREMENTS
%   none

global M_ 
for j=1:M_.endo_nbr
    assignin('caller',[deblank(M_.endo_names(j,:)) '_SS'],get_mean(deblank(M_.endo_names(j,:))));
end
