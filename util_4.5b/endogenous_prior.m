function info = endogenous_prior(T,R,Model,DynareOptions,DynareResults);

endo_priors_= DynareOptions.endo_priors_;
info=0;
infos=[0 0];
varlist=Model.endo_names(DynareResults.dr.order_var,:);
varlist=varlist(DynareResults.dr.restrict_var_list,:);
for j=1:size(endo_priors_.irf_restriction,1),
    iendo=strmatch(endo_priors_.irf_restriction{j,1},varlist,'exact');
    iexo=strmatch(endo_priors_.irf_restriction{j,2},Model.exo_names,'exact');
    if (R(iendo,iexo)>endo_priors_.irf_restriction{j,3}(1)) && (R(iendo,iexo)<endo_priors_.irf_restriction{j,3}(2)),
        infos(j,:)=[0, 0];
    else
        if R(iendo,iexo)<endo_priors_.irf_restriction{j,3}(1),
            delt = (R(iendo,iexo)-endo_priors_.irf_restriction{j,3}(1))^2;
        else
            delt = (R(iendo,iexo)-endo_priors_.irf_restriction{j,3}(2))^2;
        end            
        infos(j,:)=[201, delt];
    end
end
if any(infos),
    info=[201,sum(infos(:,2))];
end
return


