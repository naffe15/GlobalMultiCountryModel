function [rmseVECM, rmseKVECM] = vecm_prep(varargin);
global options_ bayestopt_ M_ oo_
global opts_vecm_
ys_ = oo_.steady_state;
lgy_ = M_.endo_names;
fname_ = M_.fname;

beta=opts_vecm_.beta;
gflag=opts_vecm_.gflag;

if exist(options_.datafile)
    instr = options_.datafile;
else
    instr = ['load ' options_.datafile];
end
eval(instr);
if ~exist('T'), T=[1:options_.nobs]'; end
istart = options_.presample+1;
fobs = options_.first_obs;
nobs= options_.nobs;
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
        dum=[0; diff(dum)];
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

%nex=[14 15 16];
%nk=4;

[bigVAR, yy0, yy2] = vecm0(xx, xvec, beta, opts_vecm_.nex, opts_vecm_.nk);
[bigVAR1, yy10, yy12] = vecm2(xx, xvec, beta, opts_vecm_.nex, opts_vecm_.nk, 1);
save([fname_,'_vecm'],'yy0','yy2','bigVAR','yy10','yy12','bigVAR1')

% yy2 out of sample
% yy0 within sample (plus projections)
xxx=[xx xvec];

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
        if options_.nodisplay
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
        if options_.nodisplay
            close(gcf)
        end
    end
end


% out of sample
% figure, plot([xxx(:,j)]), 
% hold on, plot([size(xx,1):size(xx,1)+nk], [xxx(end,j); yy2(:,j)],'w')
% hold on, plot([size(xx,1):size(xx,1)+nk], [xxx(end,j); yy12(:,j)],'c')
% legend('data','projection')
