ff=dir('gemc/metropolis/gemc_mh1_blck*.mat');
l1=[];
l2=[];
for k=1:length(ff)
    load([ff(k).folder filesep ff(k).name],'x2','logpo2')
    ll= nan(size(logpo2));
    for jj=1:length(logpo2)
        ll(jj) = evaluate_posterior_kernel(x2(jj,:)',M_,estim_params_,oo_,options_,bayestopt_);
        
    end
    l1=[l1 ll];
    tmp=corrcoef(ll,logpo2);
    chk(k)=tmp(2,1)-1;
    disp(chk)
    l2=[l2 logpo2];
end
