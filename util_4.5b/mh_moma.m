function th_std = mh_moma(var_list_)
global M_ oo_

DirectoryName = CheckPath('Metropolis');
load([DirectoryName,'\',M_.fname,'_param_irf1.mat'])

h = waitbar(0,'Computing MH theoretical moments...');

B=size(stock,1);
th_std=zeros(size(var_list_,1),B);
for j=1:B,
  set_all_parameters(stock(j,:))
  stoch_simul(var_list_);
  th_std(:,j)=sqrt(diag(oo_.var));
  waitbar(j/B,h)
end
close(h)


