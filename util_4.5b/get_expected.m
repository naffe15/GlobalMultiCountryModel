function y0 = get_expected(varargin);
global oo_ options_ M_

ys_ = oo_.steady_state; 
lgy_ = M_.endo_names;
ss_ = getSmootherInfo(M_,options_,oo_);
aE=ss_.aE(oo_.dr.inv_order_var,:);
% aE=ss_.aEt(oo_.dr.inv_order_var,:); use this for REAL time

lgobs_= [];
mfys=[];
for j=1:length(varargin),
    mfys = strmatch(varargin{j},lgy_,'exact');
    y0(:,j)=aE(mfys,:)'+ys_(mfys);
end

