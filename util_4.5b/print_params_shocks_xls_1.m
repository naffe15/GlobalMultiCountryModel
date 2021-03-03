function print_params_shocks_xls_1(M_,bayestopt_, CountryName)

if nargin <3
    CountryName='';
end
date_now = clock;
if date_now(5) < 10
    date_now = strcat(num2str(date_now(1)),'_',num2str(date_now(2)),'_', ...
        num2str(date_now(3)), '_at_',num2str(date_now(4)),'_', '0',num2str(date_now(5)));
    
else
    date_now = strcat(num2str(date_now(1)),'_',num2str(date_now(2)),'_', ...
        num2str(date_now(3)), '_at_',num2str(date_now(4)),'_', num2str(date_now(5)));
    
end



if isempty(bayestopt_)
    pnam=[];
else
    pnam=bayestopt_.name;
end





sd=sqrt(diag(M_.Sigma_e));

for j=1:M_.exo_nbr,
    RAW{j,2}=deblank(M_.exo_names(j,:));
    RAW{j,3}=sd(j);
    if strmatch(deblank(M_.exo_names(j,:)),pnam,'exact'),
        RAW{j,1}='**';
    else
        RAW{j,1}=' ';
    end
end


for j =1:length(M_.params)
    RAW{M_.exo_nbr+j,2}=deblank(M_.param_names(j,:));
    RAW{M_.exo_nbr+j,3}=M_.params(j);
    if strmatch(deblank(M_.param_names(j,:)),pnam,'exact'),
        RAW{M_.exo_nbr+j,1}='**';
    else
        RAW{M_.exo_nbr+j,1}=' ';
    end
end
fname=strcat(date_now,'_Parameters_Shocks.xls' );
if ~ismac
    [SUCCESS,MESSAGE]=xlswrite(fname,RAW,'Params and Shocks')
else
    [SUCCESS,MESSAGE]=xlswrite_MACOS(fname,RAW,'Params and Shocks')
end