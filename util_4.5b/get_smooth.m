
function y0 = get_smooth(varargin);
global oo_ M_

ys_ = [oo_.steady_state; zeros(M_.exo_nbr,1)];
lgy_ = char(M_.endo_names,M_.exo_names);
SmoothedVariables=[struct2cell(oo_.SmoothedVariables); struct2cell(oo_.SmoothedShocks)];
my_field_names = [fieldnames(oo_.SmoothedVariables); fieldnames(oo_.SmoothedShocks)];
isvar=zeros(length(SmoothedVariables),1);
for jf = 1:length(SmoothedVariables),
    isvar(jf)=~(isstruct(SmoothedVariables{jf}));
end
SmoothedVariables=cell2struct(SmoothedVariables(logical(isvar)),my_field_names(logical(isvar)));


lgobs_= [];
mfys=[];
y0=zeros(length(eval(['SmoothedVariables.',varargin{1}])),length(varargin));
for j=1:length(varargin),
%     mfys = strmatch(varargin{j},lgy_,'exact');
    y0(:,j)=eval(['SmoothedVariables.',varargin{j}]); %+ys_(mfys);
end

