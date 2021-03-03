function append_mode_values(fname, pname, pval)
% append_mode_values(fname, pname, pval)
% inputs
% fname = char name of the mode file
% pname = cell array of estimated param of shock
% pval = real array of new values for mode files

if nargin==0,
    disp('append_mode_values(fname, pname, pval)')
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
    indx=strmatch(pname{j},parameter_names,'exact');
    if isempty(indx),
        parameter_names = [parameter_names; pname(j)];
        xparam1 = [xparam1; pval(j)];
    else
        error([pname{j} 'is already in mode file!'])
    end
end
hh=eye(length(xparam1));
save(fname,'hh','xparam1','parameter_names','-append')

end

