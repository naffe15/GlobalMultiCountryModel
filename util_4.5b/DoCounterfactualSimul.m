function [my_simul , my_newalpha] = DoCounterfactualSimul(myetahat,myalphahat,t0,oo,M,options,initial_condition)
% DoCounterfactutalSimul This function runs a counterfactual simulation

% INPUTS:   myetahat            = (nexo,periods) matrix of shocks 
%           initial_condition   = (switch for initial conditions)
% OUTPUTS:  simul_mat           = structure of endogeneous vars (#nendo)
%                                 each for the smoothed time periods
if nargin < 4
    global M_ oo_ options_
    M = M_;
    oo = oo_;
    options = options_;
else
    
end
if nargin < 7
    initial_condition = 1;
end

if isempty(myalphahat)
	ss_ = getSmootherInfo(M,options,oo);
	myalphahat = ss_.a1;
end
%ss_ = getSmootherInfo(M,options,oo);
% oo  = load('smoother_info', 'oo_'); % in case ordering of variables is
% different.

%% Resolve CC model to retrieve A, B, and so on
[A,B,ys,info] = dynare_resolve(M, options, oo);

%% Counterfactual simulations with initial condition
if initial_condition
    y0(oo.dr.inv_order_var)= myalphahat(oo.dr.inv_order_var,t0);
    %ss_.a1(oo_.dr.inv_order_var);
    % y0(oo_.dr.inv_order_var)=ss_.a1(oo.oo_.dr.inv_order_var);
    y2=y0';
else
    ss_ = getSmootherInfo(M,options,oo);
    y2=ss_.a1*0;
end

shocks = myetahat(:,t0+1:end);

for j=1:size(shocks,2)
    y2(:,j+1)=A*y2(:,j)+B*shocks(:,j);
end

% Re-ordering ofvariables 
my_newalpha =y2(oo.dr.inv_order_var,:);


%% Store simulation

my_simul = struct();
for i = 1:M.endo_nbr
    my_var = deblank(M.endo_names(i,:));
    my_ind = varlist_indices(my_var,M.endo_names);
    my_simul.(my_var) = my_newalpha(my_ind,:);
end

% Re-ordering ofvariables 
my_newalpha =my_newalpha(oo.dr.order_var,:);

disp('Counterfactural simul completed')
end

