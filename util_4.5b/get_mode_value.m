function x = get_mode_value(fname,pname)
% get_mode_values(fname, pname)
% inputs
% fname = char name of the mode file
% pname = cell array of estimated param of shock

global M_
icheck=load(fname,'parameter_names');
if isempty(icheck),
    disp('parameter_names is not present!')
else
    parameter_names = icheck.parameter_names;
end
load(fname,'xparam1')

i = strmatch(pname,parameter_names,'exact');

if isempty(i)
    error(sprintf('Can''t find parameter %s', pname))
end

x = xparam1(i);
