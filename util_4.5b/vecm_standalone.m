function [out,rmseVECM, rmseKVECM] = vecm_standalone(opts_,varargin)
global opts_vecm_

opts_vecm_ = opts_;

fname_ = opts_vecm_.fname_;
beta=opts_vecm_.beta;
gflag=opts_vecm_.gflag;
nlags=opts_vecm_.nlags;

if exist(opts_vecm_.datafile)
    instr = opts_vecm_.datafile;
else
    instr = ['load ' opts_vecm_.datafile];
end
eval(instr);
if ~exist('T'), T=[1:opts_vecm_.nobs]'; end
istart = opts_vecm_.presample+1;
fobs = opts_vecm_.first_obs;
nobs= opts_vecm_.nobs;
if any(gflag) & opts_vecm_.lflag,
    istart=istart-1;
    fobs=fobs+1;
    nobs=nobs-1;
end

xx=[];
x0=[];
for j=1:length(varargin),
    eval(['dum=',varargin{j},';'])
    x0=[x0 dum];
    if gflag(j),
        dum=[NaN; diff(dum)];
    end
    xx=[xx dum];
end
%xx = [E_GC  E_GE  E_GEX  E_GI  E_GIM  E_GL  E_GY  E_WPHI  E_INOM  ...
%     E_PHI  E_PHIC  E_PHIM  E_PHIX  E_GTR  E_INOMW  E_PHIW  E_GYW];  
xx=xx(fobs:fobs+nobs-1,:);
     
xvec=[];
for j=1:size(beta,1),
    xvec = [xvec x0*beta(j,:)'];
end
%xvec = [ E_LCSN  E_LISN  E_LEXY  E_LIMY  E_LPCP  E_LPMP  E_LPXP  E_L  E_LWRY  E_LTRYN  E_LYWY];
if ~isempty(xvec)
    xvec=xvec(fobs:fobs+nobs-1,:);
end
% check cointegrations;
figure('name','Check cointegration plots'),
for j=1:size(beta,1),
    element = find(beta(j,:));
    tit=[];
    for ielem=1:length(element);
        tit=[tit, int2str(beta(j,element(ielem))),'*',varargin{element(ielem)}(3:end),'  '];
    end
    subplot(4,4,j), plot(xvec(:,j)), title(tit,'interpreter','none'),
end


%nex=[14 15 16];
%nk=4;

% [bigVAR, yy0, yy2, ym] = vecm0(xx, xvec, beta, opts_vecm_.nex, opts_vecm_.nk);
[bigVAR , yy0, yy2, ym] = vecm1(xx, xvec, beta, opts_vecm_.nex, opts_vecm_.nk, nlags);
[bigVAR1, yy10, yy12, ym1] = vecm2(xx, xvec, beta, opts_vecm_.nex, opts_vecm_.nk, nlags);

% yy2 out of sample
% yy0 within sample (plus projections)
xxx=[xx xvec];
save([fname_,'_vecm'],'xxx','yy0','yy2','bigVAR','yy10','yy12','bigVAR1')

ifig=0;disp(' ')
%disp(' ')
%disp('            VECM RMSE')
for j=1:length(varargin),
    rmseVECM(j)=sqrt(mean((xxx(istart:end,j)-squeeze(yy0(istart:size(xx,1),j,1))).^2)); 
    rmseVECM1(j)=sqrt(mean((xxx(istart:end,j)-squeeze(yy10(istart:size(xx,1),j,1))).^2)); 
    if gflag(j)
        ttit=['diff(',varargin{j},')'];
    else
        ttit=varargin{j};
    end
    %disp([ttit, sprintf('%15.5g',[rmseVECM(j)'])])
    if mod(j,9)==1,
        figure('Name',['VECM 1-step ahead prediction']),
        ifig=ifig+1;
    end
    subplot(3,3,j-9*(ifig-1))
    plot(T(fobs+istart-1:fobs+nobs-1),[xxx(istart:end,j) squeeze(yy0(istart:size(xx,1),j,1)) squeeze(yy10(istart:size(xx,1),j,1))]) % example within sample plot
    title(ttit,'interpreter','none')
    if mod(j,9)==0 | j==length(varargin),
        saveas(gcf,[fname_,'_VECM_Fit',int2str(ifig)])
        if opts_vecm_.nodisplay
            close(gcf)
        end
    end
end

nk=opts_vecm_.nk;
ifig=0;
%disp(' ')
%disp(' ')
%disp('            VECM k-step RMSE')
for j=1:length(varargin),
    rmseKVECM(j)=sqrt(mean((xxx(nk+1:end,j)-squeeze(yy0(nk+1:size(xx,1),j,nk))).^2)); 
    rmseKVECM1(j)=sqrt(mean((xxx(nk+1:end,j)-squeeze(yy10(nk+1:size(xx,1),j,nk))).^2)); 
    if gflag(j)
        ttit=['diff(',varargin{j},')'];
    else
        ttit=varargin{j};
    end
    %disp([ttit, sprintf('%15.5g',[rmseKVECM(j)'])])
    if mod(j,9)==1,
        figure('Name',['VECM ',int2str(nk),'-step ahead prediction']),
        ifig=ifig+1;
    end
    subplot(3,3,j-9*(ifig-1))
    plot(T(fobs+nk:fobs+nobs-1),[xxx(nk+1:end,j) squeeze(yy0(nk+1:size(xx,1),j,nk)) squeeze(yy10(nk+1:size(xx,1),j,nk))]) % example within sample plot
    title(ttit,'interpreter','none')
    if mod(j,9)==0 | j==length(varargin),
        saveas(gcf,[fname_,'_VECM_',num2str(nk),'-step',int2str(ifig)])
        if opts_vecm_.nodisplay
            close(gcf)
        end
    end
end

out.xxx=xxx;
out.bigVAR=bigVAR;
out.bigVAR1=bigVAR1;
out.yy0=yy0;
out.yy2=yy2;
out.yy10=yy10;
out.yy12=yy12;
out.rmseVECM=rmseVECM;
out.rmseVECM1=rmseVECM1;
out.rmseKVECM=rmseKVECM;
out.rmseKVECM1=rmseKVECM1;

ifig=0;
Tend = T(fobs+nobs-1)+opts_vecm_.forecast*0.25;
for j=1:length(varargin),
    if gflag(j)
        ttit=['diff(',varargin{j},')'];
    else
        ttit=varargin{j};
    end
    if mod(j,9)==1,
        figure('Name',['VECM projections']),
        ifig=ifig+1;
    end
    subplot(3,3,j-9*(ifig-1))
    plot(T(fobs+nobs-80:fobs+nobs-1),xxx(end-79:end,j)),
    hold on,
    plot(T(fobs+nobs-1)+0.25:0.25:Tend,[yy12(:,j)],'.-r') 
    plot(T(fobs+nobs-1)+0.25:0.25:Tend,[yy2(:,j)],'.-k') 
    plot([T(fobs+nobs-80),Tend],[ym1(j) ym1(j)],':r') 
    plot([T(fobs+nobs-80),Tend],[ym(j) ym(j)],':k') 
%     legend('data','VECM','VECM with theo mean','VECM mean','vecm THEO mean')
    title(ttit,'interpreter','none')
    if mod(j,9)==0 | j==length(varargin),
        saveas(gcf,[fname_,'_VECM_projections',int2str(ifig)])
        if opts_vecm_.nodisplay,
            close(gcf)
        end
    end
end



