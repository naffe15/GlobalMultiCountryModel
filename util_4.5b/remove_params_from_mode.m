function remove_params_from_mode(fname, pname)
% remove_params_from_mode(fname, pname)
% inputs
% fname = char name of the mode file
% pname = cell array of estimated param of shock

if nargin==0,
    disp('remove_params_from_mode(fname, pname)')
    return
end

icheck=load(fname,'parameter_names');
if isempty(icheck),
    disp('parameter_names is not present!')
else
    parameter_names = icheck.parameter_names;
end
load(fname,'xparam1','hh')
for j=1:length(pname)
    indx(j)=strmatch(pname{j},parameter_names);
end
xparam1(indx)=[];
parameter_names(indx)=[];
hh(indx,:)=[];
hh(:,indx)=[];

save(fname,'xparam1','hh','parameter_names')

end

