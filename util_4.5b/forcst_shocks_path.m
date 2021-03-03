% function [yf,int_width]=forcst_shocks_path(dr,y0,horizon,var_list,shocks_path)
function info = forcst_shocks_path(dr,y0,horizon,var_list,shocks_path)

global M_  oo_ options_ 

info = 0;

if nargin<5 || isempty(shocks_path)
    shocks_path=zeros(horizon,M_.exo_nbr);
end


make_ex_;


maximum_lag = M_.maximum_lag;

endo_names = M_.endo_names;
if isempty(var_list)
    var_list = endo_names(1:M_.orig_endo_nbr, :);
end
i_var = [];
for i = 1:size(var_list)
    tmp = strmatch(var_list(i,:),endo_names,'exact');
    if isempty(tmp)
        error([var_list(i,:) ' isn''t and endogenous variable'])
    end
    i_var = [i_var; tmp];
end

n_var = length(i_var);



yf = simult_(y0,dr,shocks_path,1);
nstatic = M_.nstatic;
nspred = M_.nspred;
nc = size(dr.ghx,2);
endo_nbr = M_.endo_nbr;
inv_order_var = dr.inv_order_var;
[A,B] = kalman_transition_matrix(dr,nstatic+(1:nspred),1:nc,M_.exo_nbr);

if size(var_list,1) == 0
    var_list = M_.endo_names(1:M_.orig_endo_nbr,:);
end
nvar = size(var_list,1);
ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
    if isempty(i_tmp)
        disp(var_list(i,:));
        error (['One of the variable specified does not exist']) ;
    else
        ivar(i) = i_tmp;
    end
end

ghx1 = dr.ghx(inv_order_var(ivar),:);
ghu1 = dr.ghu(inv_order_var(ivar),:);

sigma_u = B*M_.Sigma_e*B';
sigma_u1 = ghu1*M_.Sigma_e*ghu1';
sigma_y = 0;

var_yf=NaN(horizon,nvar); %initialize
for i=1:horizon
    sigma_y1 = ghx1*sigma_y*ghx1'+sigma_u1;
    var_yf(i,:) = diag(sigma_y1)';
    if i == horizon
        break
    end
    sigma_u = A*sigma_u*A';
    sigma_y = sigma_y+sigma_u;
end

fact = norminv((1-options_.conf_sig)/2,0,1);

int_width = zeros(horizon,M_.endo_nbr);
for i=1:nvar
    int_width(:,i) = -fact*sqrt(var_yf(:,i));
end

yf = yf(ivar,:);



for i=1:n_var
    eval(['oo_.forecast_exo_path.Mean.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+(1:horizon))'';']);
    eval(['oo_.forecast_exo_path.HPDinf.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+(1:horizon))''-' ...
          ' int_width(1:horizon,' int2str(i) ');']);
    eval(['oo_.forecast_exo_path.HPDsup.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+(1:horizon))''+' ...
          ' int_width(1:horizon,' int2str(i) ');']);
end


