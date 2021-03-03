function [ax, etahat] = histo_simul(T, R, a0, etahat, exo_info),
global M_

if nargin<5, exo_info=[]; end
% exo_flag=ones(M_.exo_nbr,1);
if ~isempty(exo_info),
  for j=1:size(exo_info.names,1),
    iexo = strmatch(exo_info.names(j,:),M_.exo_names,'exact');
    etahat(iexo,exo_info.t0:end)=exo_info.value(j);
  end
end
  

ax(:,1) = a0;
for j=1:size(etahat,2), 
  ax(:,j+1)=T*ax(:,j)+R*etahat(:,j); 
end

