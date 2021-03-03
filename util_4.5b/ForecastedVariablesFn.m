function out=ForecastedVariablesFn(nfrcst, M_,options_,oo_, addsteady, exnam, exval)



outCC=CrossCovariance(M_, oo_);
ccova=outCC.ccova;

ss_ = getSmootherInfo(M_, options_, oo_);

ex=zeros(size(ss_.R,2), nfrcst);
if nargin==7,
    for j=1:length(exnam),
        ix=strcmp(exnam{j},cellstr(M_.exo_names));
        ex(ix,:)=exval(:,j);
    end
end

dsteady = oo_.steady_state-ss_.SteadyState;
yfcst=ss_.a+dsteady(oo_.dr.order_var);
Vfrcst=zeros(length(yfcst), length(yfcst), nfrcst);
for j=1:nfrcst,
    yfcst(:,j+1)=ss_.T*yfcst(:,j)+ss_.R*ex(:,j);
    Vfrcst(:,:,j+1)=ss_.T*Vfrcst(:,:,j)*(ss_.T)'+ss_.R*ccova*(ss_.R)';
end
fnam = fieldnames(oo_.SmoothedVariables);


if addsteady
    for j=1:length(fnam),
        s(j)=get_mean(fnam{j});
        eval(['ForecastedVariables.',fnam{j},'= s(j)+[oo_.SmoothedVariables.',fnam{j},'(1:end)-s(j)+dsteady(oo_.dr.order_var(j)); yfcst(j,2:end)''];'])
%         eval(['ForecastedVariables.',fnam{j},'=get_mean(',fnam{j},')+[oo_.SmoothedVariables.',fnam{j},'(1:end); yfcst(j,2:end)''];'])
        eval(['ForecastedVariablesStd.',fnam{j},'=[oo_.SmoothedVariables.',fnam{j},'(1:end)*NaN; squeeze(sqrt(Vfrcst(j,j,2:end)))];'])
    end
else
    for j=1:length(fnam),    
        s(j)=get_mean(fnam{j});
        eval(['ForecastedVariables.',fnam{j},'=[oo_.SmoothedVariables.',fnam{j},'(1:end)+dsteady(oo_.dr.order_var(j))-s(j); yfcst(j,2:end)''];'])
        eval(['ForecastedVariablesStd.',fnam{j},'=[oo_.SmoothedVariables.',fnam{j},'(1:end)*NaN; squeeze(sqrt(Vfrcst(j,j,2:end)))];'])
        
    end
    
end
out.ForecastedVariables=ForecastedVariables;
out.ForecastedVariablesStd=ForecastedVariablesStd;
out.s=s;

save([M_.fname,'_frcst'],'out');