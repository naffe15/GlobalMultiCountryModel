function y0 = getFilterRuns(varlist)
global M_ oo_

y0=[];
dirname = CheckPath('Metropolis');
filfilt = dir([dirname '/' M_.fname '_filter*.mat']);
order_var=oo_.dr.order_var;
for j=1:size(varlist,1)
   jxj(j)=strmatch(deblank(varlist(j,:)),M_.endo_names,'exact');
   jfilt(j)=find(order_var==jxj(j));
end

filparam = dir([dirname '/' M_.fname '_param*.mat']);
sto_ys=[];
for j=1:length(filparam),
  %load([DirectoryName '/' M_.fname '_param',int2str(j),'.mat']);
  if isempty(strmatch([M_.fname '_param_irf'],filparam(j).name))
    load([dirname '/' filparam(j).name],'stock_ys');
    sto_ys=[sto_ys; stock_ys];
    clear stock stock_logpo stock_ys;
  end
end


nbb=0;
for j=1:length(filfilt),
  load([dirname '/' M_.fname '_filter',num2str(j),'.mat']);
  nb = size(stock,4);
  %y0(:,nbb+1:nbb+nb)=squeeze(stock(1,jxj,:,:));
  for i=1:length(jxj)
%     y0(i,:,nbb+1:nbb+nb)=squeeze(stock(1,jfilt(i),:,:)) + ...
    y0(i,:,nbb+1:nbb+nb)=squeeze(stock(1,jxj(i),:,:)) + ...
      kron(sto_ys(nbb+1:nbb+nb,jxj(i))',ones(size(stock,3),1));
  end
  %y0(:,:,size(y0,3):size(y0,3)+size(stock,3))=stock;
  nbb=nbb+nb;
end
