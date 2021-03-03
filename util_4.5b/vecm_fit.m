function [rmse, r2, vecm]=vecm_fit(nlags, iplot),
% [rmse, r2, vecm]=vecm_fit(nlags, iplot),

global options_ opts_vecm_ M_

fname_ = M_.fname;

if nargin==0,
  disp(' ')
  disp('[rmse, r2, vecm]=vecm_fit(nlags, iplot)')
  return
end

if nargin<2,
  iplot=0;
end

[bvar, nyy, nxx] = bvar_fit(nlags);

bigvar=bvar.B';
nstate=size(bigvar,2);
for j=1:nlags-1,
  dum=zeros(nyy,nstate);
  dum(:,[1:nyy]+nyy*(j-1))=eye(nyy);
  bigvar=[bigvar; dum];
end
ee=eig(bigvar(:,1:end-1));
if max(abs(ee))>1,
  warning('The estimated BVAR is unstable!')
end

[vecm, ny, nx, forecast_data] = vecm_toolbox(nlags);

if ~isempty(vecm)
bigrho=vecm.B';
beta=vecm.beta;
nstate=size(bigrho,2);

dumVECM=[];
for j=1:nlags,
    dumVECM = [dumVECM beta*bigrho(:,1+ny*(j-1):ny*j)];
end
dumVECM=[dumVECM  eye(nx)+beta*bigrho(:,[ny*nlags+1:end]) ];

for j=1:nlags-1,
  dum=zeros(ny,nstate);
  dum(:,[1:ny]+ny*(j-1))=eye(ny);
  bigrho=[bigrho; dum];
end

bigrho=[bigrho; dumVECM];

ee=eig(bigrho);
if max(abs(ee))>1,
  warning('The estimated VECM is unstable!')
end
vecm.BB=bigrho;

else
  vecm=bvar;
  vecm.BB=bigvar;
  idx = options_.first_obs+options_.presample:options_.first_obs+options_.nobs-1;
  no=length(idx);
  vecm.u=vecm.u(1:no,:,:);
  vecm.yhat=vecm.yhat(1:no,:,:);
end
  
for j=1:size(vecm.u,3),
  rmse(:,j)=sqrt(mean(squeeze(vecm.u(j:end,:,j)).^2));
  r2(:,j)=1-var(squeeze(vecm.u(j:end,:,j)))./var(vecm.ydata(nlags+j:end,:));
end

disp('1-step  ahead predictions of VECM')
disp('             RMSE           R2')
% vname=str2mat(varargin{:});
for j=1:ny,
%     iv = strmatch(deblank(lgobs_(j,:)), lgobs_,'exact');
    disp([opts_vecm_.varobs(j,:), sprintf('%15.5g',[rmse(j,1)' r2(j,1)])])
end,

nk=size(vecm.u,3);
disp(' ')
disp([int2str(nk),'-step  ahead predictions of VECM'])
disp('             RMSE           R2')
% vname=str2mat(varargin{:});
for j=1:ny,
%     iv = strmatch(deblank(lgobs_(j,:)), lgobs_,'exact');
    disp([opts_vecm_.varobs(j,:), sprintf('%15.5g',[rmse(j,end)' r2(j,end)])])
end,

if iplot,
ifig=0;
for j=1:ny,
    if mod(j,9)==1,
        figure('Name',['VECM 1-step ahead prediction']),
        ifig=ifig+1;
    end
    subplot(3,3,j-9*(ifig-1))
    plot(vecm.T(nlags+1:end),[squeeze(vecm.yhat(:,j,1)) vecm.ydata(nlags+1:end,j)]) % example within sample plot
    hh=get(gca,'children');
    set(hh(end-1),'linestyle','-','color',[0.7 0.7 0.7])
    set(hh(end),'linestyle','-','color',[0 0 0])
    title(opts_vecm_.varobs(j,:),'interpreter','none')
    if mod(j,9)==0 | j==ny,
        saveas(gcf,[fname_,'_VECM_Fit',int2str(ifig)])
        if options_.nodisplay
            close(gcf)
        end
    end
end

if nk>1,
ifig=0;
for j=1:ny,
    if mod(j,9)==1,
        figure('Name',['VECM ',int2str(nk),'-step ahead prediction']),
        ifig=ifig+1;
    end
    subplot(3,3,j-9*(ifig-1))
    plot(vecm.T(nlags+nk:end),[squeeze(vecm.yhat(nk:end,j,nk)) vecm.ydata(nlags+nk:end,j)]) % example within sample plot
    title(opts_vecm_.varobs(j,:),'interpreter','none')
    if mod(j,9)==0 | j==ny,
        saveas(gcf,[fname_,'_VECM_Fit_K',int2str(ifig)])
        if options_.nodisplay
            close(gcf)
        end
    end
end
end
  
  
end




