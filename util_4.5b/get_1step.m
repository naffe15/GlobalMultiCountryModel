function y0 = get_1step(varargin);
global oo_ M_

ys_ = oo_.steady_state; 
lgy_ = M_.endo_names;
ss_ = getSmootherInfo(oo_);
alphat=ss_.alphat(oo_.dr.inv_order_var,:);

lgobs_= [];
mfys=[];
for j=1:length(varargin),
    mfys = strmatch(varargin{j},lgy_,'exact');
    y0(:,j)=alphat(mfys,:)'+ys_(mfys);
end

