function [my_simul] = DoCounterfactutalSimul(myetahat,initial_condition)
% DoCounterfactutalSimul This function runs a counterfactual simulation

% INPUTS:   myetahat            = (nexo,periods) matrix of shocks 
%           initial_condition   = (switch for initial conditions)
% OUTPUTS:  simul_mat           = structure of endogeneous vars (#nendo)
%                                 each for the smoothed time periods

if nargin < 2
    initial_condition = 1;
end
    
global M_ oo_ options_

ss_ = getSmootherInfo(M_,options_,oo_);
% oo  = load('smoother_info', 'oo_'); % in case ordering of variables is
% different.

%% Resolve CC model to retrieve A, B, and so on
[A,B,ys,info] = dynare_resolve(M_, options_, oo_);

%% Counterfactual simulations with initial condition
if initial_condition
    y0(oo_.dr.inv_order_var)=ss_.a1(oo_.dr.inv_order_var);
    % y0(oo_.dr.inv_order_var)=ss_.a1(oo.oo_.dr.inv_order_var);
    y2=y0';
else
    y2=ss_.a1*0;
end
    
for j=1:size(ss_.etahat,2)-1
    y2(:,j+1)=A*y2(:,j)+B*myetahat(:,j+1);
end
% Re-ordering of variables 
simul_mat=y2(oo_.dr.inv_order_var,:);

%% Store simulation
my_simul = struct();
for i = 1:M_.endo_nbr
    my_var = deblank(M_.endo_names(i,:));
    my_ind = varlist_indices(my_var,M_.endo_names);
    my_simul.(my_var) = simul_mat(my_ind,:);
end

disp('Counterfactural simul completed')
end

