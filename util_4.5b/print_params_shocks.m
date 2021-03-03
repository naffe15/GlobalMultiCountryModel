function print_params_shocks;

global M_ bayestopt_

if isempty(bayestopt_)
  pnam=[];
else
  if isfield(bayestopt_,'name')
      pnam=bayestopt_.name;
  else
      pnam=[];
  end
end
  

disp('PARAMETER VALUES (** for estimated parameters)')
for j=1:length(M_.params)
  if strmatch(deblank(M_.param_names(j,:)),pnam,'exact'),
    disp(['**  ',M_.param_names(j,:),'  ',num2str(M_.params(j),'%21.14g')])
    
  else
    
    disp(['    ',M_.param_names(j,:),'  ',num2str(M_.params(j),'%21.14g')])
  end
end

disp(' ')
disp('SHOCKS STANDARD DEVIATION (** for estimated shocks)')
sd=sqrt(diag(M_.Sigma_e));
for j=1:M_.exo_nbr,
  if strmatch(deblank(M_.exo_names(j,:)),pnam,'exact'),
    disp(['**  ',M_.exo_names(j,:),'  ',num2str(sd(j),'%21.14g')])
    
  else
    
    disp(['    ',M_.exo_names(j,:),'  ',num2str(sd(j),'%21.14g')])
  end
end