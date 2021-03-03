function set_mode_values(fname, pname, pval)
% set_mode_values(fname, pname, pval)
% inputs
% fname = char name of the mode file
% pname = cell array of estimated param of shock
% pval = real array of new values for mode files

if nargin==0,
    disp('set_mode_values(fname, pname, pval)')
    return
end

icheck=load(fname,'parameter_names');
if isempty(icheck),
    disp('parameter_names is not present!')
else
    parameter_names = icheck.parameter_names;
end
load(fname,'xparam1')
for j=1:length(pname)
    indx=strmatch(pname{j},parameter_names);
    xparam1(indx) = pval(j);
end
save(fname,'xparam1','-append')

end

