function [VVi, VVinf, th_std] = mh_vardec(k,varargin)
global M_ oo_

DirectoryName = CheckPath('Metropolis');
load([DirectoryName,'\',M_.fname,'_param_irf1.mat'])

h = waitbar(0,'Computing MH variance decomposition ...');
var_list_=str2mat(varargin);

B=size(stock,1);
th_std=zeros(size(var_list_,1),B);
VVi=zeros(size(var_list_,1),M_.exo_nbr,k,B);
VVinf=zeros(size(var_list_,1),M_.exo_nbr,B);
for j=1:B,
  set_all_parameters(stock(j,:))
  stoch_simul(var_list_);
  th_std(:,j)=sqrt(diag(oo_.var));
  [Vi, Vinf]=vardec(k,varargin{:});
  VVi(:,:,:,j)=Vi;
  VVinf(:,:,j)=Vinf;

  waitbar(j/B,h)
end
close(h)


