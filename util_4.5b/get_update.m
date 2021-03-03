function y0 = get_update(varargin);
global oo_ M_

ys_ = oo_.steady_state; 
lgy_ = M_.endo_names;


lgobs_= [];
mfys=[];
y0=zeros(length(eval(['oo_.UpdatedVariables.',varargin{1}])),length(varargin));
for j=1:length(varargin),
%     mfys = strmatch(varargin{j},lgy_,'exact');
    y0(:,j)=eval(['oo_.UpdatedVariables.',varargin{j}]); %+ys_(mfys);
end

