function out=SimulatedVariablesFn(nfrcst, M_,options_,oo_, addsteady, exnam, exval)


[T,R,SteadyState] = dynare_resolve(M_, options_, oo_);

ex=zeros(size(R,2), nfrcst);
if nargin==7,
    for j=1:length(exnam),
        ix=strcmp(exnam{j},cellstr(M_.exo_names));
        ex(ix,:)=exval(:,j);
    end
end

dsteady = oo_.steady_state-SteadyState;
yfcst=zeros(M_.endo_nbr,1)+dsteady(oo_.dr.order_var);
for j=1:nfrcst,
    yfcst(:,j+1)=T*yfcst(:,j)+R*ex(:,j);
end
fnam = cellstr(M_.endo_names(oo_.dr.order_var,:));


if addsteady
    for j=1:length(fnam),
        s(j)=get_mean(fnam{j});
        eval(['SimulatedVariables.',fnam{j},'= s(j)+[yfcst(j,2:end)''];'])
    end
else
    for j=1:length(fnam),    
        s(j)=get_mean(fnam{j});
        eval(['SimulatedVariables.',fnam{j},'=[yfcst(j,2:end)''];'])
       
    end
    
end
out.SimulatedVariables=SimulatedVariables;
out.s=s;

save([M_.fname,'_simul'],'out');