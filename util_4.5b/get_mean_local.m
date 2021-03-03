function y0 = get_mean_local(M_,oo_,options_,varargin);
 

if ~isempty(regexp(varargin{end},'\d')) && isempty(regexp(varargin{end},'\D')),
    order=eval(varargin{end});
else
    order=1;
end
if order==1,
     ys_ = oo_.steady_state;
% % %      if isfield(oo_.dr,'ys')
% % %          if length(ys_) == length(oo_.dr.ys),
% % %              ys_ = oo_.dr.ys;       %changed back to previous line due to issues
% % %                                      wrt get_mean of LKTOT  %hohbest
% % %          end
% % %      end
     [ys_,params,info] = evaluate_steady_state(ys_,M_,options_,oo_,1);
elseif order==2,
    ys_ = oo_.dr.ys;
    ys_(oo_.dr.order_var)=ys_(oo_.dr.order_var)+oo_.dr.ghs2./2;
else
   return
end
lgy_ = M_.endo_names;


lgobs_= [];
mfys=[];
for j=1:length(varargin),
    dum = strmatch(varargin{j},lgy_,'exact');
    mfys = [mfys dum];
end

y0 = ys_(mfys);
