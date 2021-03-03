fname = 'quest3hlmr';
firstchain=1;
nchain = 4;
aa=dir([fname,'/metropolis/',fname,'_mh*_blck*.mat']);
firstifile = 1;
nfiles = floor(length(aa)/nchain);
firstifile = max(1,floor(nfiles/2));
ll=[];
for j=firstchain:firstchain+nchain-1, 
    for i=firstifile:nfiles, 
    load([fname,'/metropolis/',fname,'_mh',int2str(i),'_blck',int2str(j)]); ll=[ll; logpo2]; 
end, end
nr=size(ll,1)/nchain;
figure('position',[190   160   560   750])
for j=1:length(bayestopt_.name), 
    xx=[]; 
    for jj=firstchain:firstchain+nchain-1, 
        for i=firstifile:nfiles, 
            load([fname,'/metropolis/',fname,'_mh',int2str(i),'_blck',int2str(jj)]); 
            xx=[xx; x2(:,j)]; 
        end, 
    end,
    [N,X] = hist(xx,30); [d,m]=max(N); MCMCmode(j,1)=X(m);
    subplot(311), 
    plot(reshape(xx,[nr,nchain])), 
    hold on, plot([0 nr],xparam1([j j]),'k-'), hold off,
    set(gca,'ylim',[min(xx) max(xx)]),
    
    subplot(312), 
    hist(xx,30), 
    hold on, plot(xparam1([j j]),get(gca,'ylim'),'r'), 
    hold off, title(bayestopt_.name{j}), 
    
    subplot(313), 
    plot(xx,ll,'.'), 
    [lm, jm]=max(ll);
    hold on, plot(xx(jm), ll(jm), 'r*'),
    hold off,
    set(gca,'ylim',[min(ll) max(ll)]), 
    pause; %(0.5), 
end,

figure, 
subplot(211), bar([abs((MCMCmode-xparam1)./xparam1 )])
subplot(212), bar(MCMCmode./xparam1)
xp0=xparam1;
xparam1=MCMCmode;
hh=eye(length(bayestopt_.name));
save([fname,'_mode_hist'],'xparam1','hh')
xparam1=xp0;

