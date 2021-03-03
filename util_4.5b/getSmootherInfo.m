function [ss_, DynareOutput] = getSmootherInfo(M_, options_, oo_);


[T,R,SteadyState,info,Model,DynareOptions,DynareOutput] = dynare_resolve(M_, options_, oo_);
SS = SteadyState(oo_.dr.order_var);
fnam = fieldnames(oo_.SmoothedVariables);
SmoothedVariables=[struct2cell(oo_.SmoothedVariables)];
isvar=zeros(length(SmoothedVariables),1);
for jf = 1:length(SmoothedVariables),
    isvar(jf)=~(isstruct(SmoothedVariables{jf}));
end
fnam=fnam(logical(isvar));

for j=1:length(T),
    alphahat(j,:)=getfield(oo_.SmoothedVariables,fnam{j})-SS(j);
end
for j=1:length(T),
    alphatt(j,:)=getfield(oo_.UpdatedVariables,fnam{j})-SS(j);
end

if isfield(oo_,'FilteredVariables')
    
for j=1:length(T),
    alphat(j,:)=getfield(oo_.FilteredVariables,fnam{j})-SS(j);
end

end
fnam = fieldnames(oo_.SmoothedShocks);
SmoothedShocks=[struct2cell(oo_.SmoothedShocks)];
isvar=zeros(length(SmoothedShocks),1);
for jf = 1:length(SmoothedShocks),
    isvar(jf)=~(isstruct(SmoothedShocks{jf}));
end
fnam=fnam(logical(isvar));
for j=1:size(R,2),
    etahat(j,:)=getfield(oo_.SmoothedShocks,fnam{j});
    if isfield(oo_, 'UpdatedShocks')
        eta(j,:)=getfield(oo_.UpdatedShocks,fnam{j});
    end
end

aEt=T*alphatt;
aE=T*alphahat;
ss_.T=T;
ss_.R=R;
ss_.SteadyState=SteadyState;

ss_.a=alphahat(:,end);
ss_.a1=alphahat(:,1);
ss_.ss=SteadyState(oo_.dr.order_var);
for j= (length(ss_.ss)+1):length(T),
    ss_.ss(j)=ss_.ss(find(T(j,:)));
end
ss_.etahat=etahat;
ss_.alphahat=alphahat;
ss_.alphatt=alphatt;
% ss_.af=af;
ss_.aE=aE;
ss_.aEt=aEt;

if isfield(oo_,'FilteredVariables')

ss_.alphat=alphat(:,1:size(alphahat,2));

end
if isfield(oo_, 'UpdatedShocks')
    ss_.eta=eta;
end
