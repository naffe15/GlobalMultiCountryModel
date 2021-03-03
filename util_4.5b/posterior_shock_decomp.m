function posterior_shock_decomp(varargin)
global M_ oo_ options_

load data T;
TT=[T(options_.first_obs:end); 2009];

vname = {'E_GY', 'E_GI', 'E_GIHOUSE', 'E_GIHOUSECC','E_GIHOUSENLC', 'E_INFHOUSE', ...
    'E_GC', 'E_GCCC', 'E_GCNLC', 'E_GEX', 'E_GIM', 'E_GDEBTCC','E_INFLAND','E_INFCONSTR','E_INFC'};

ex_names_={'E_EPS_LTFP',  'E_UPI';
  'E_EPS_M','';
  'E_EPS_RPREMK','';
  'E_EPS_RPREMHOUSECC',   'E_EPS_RPREMLANDE';
  'E_EPS_DEBTCCT',''};

% ex_names_={'E_EPS_LTFP';
%   'E_UPI';
%   'E_EPS_M';
%   'E_EPS_RPREME';
%   'E_EPS_RPREMK';
%   'E_EPS_RPREMHOUSECC';
%   'E_EPS_RPREMLANDE';
%   'E_EPS_XI'};
% 
% leg=[ex_names_;{'others'}];
leg={'Technology'; 'Monetary policy shock'; 'Stock market bubble'; 'Housing bubble'; 'Collateral shock'; 'Others'};

t1decomp =22; % set it manually to change the initial point for plotting shock decomp
kpoint1 = t1decomp:options_.nobs-3;  
% kpoint1 = 62:101;
kpoint2 = kpoint1+1;
kpoint3 = kpoint1+2;
kpoint4 = kpoint1+3;


varlist = char(varargin);
DirectoryName = CheckPath('metropolis');

load([DirectoryName '/' M_.fname '_data.mat'],'stock_gend','stock_data');

ffil = dir([DirectoryName '/' M_.fname '_smooth*.mat']);
fileParam = dir([DirectoryName '/' M_.fname '_param*.mat']);
fileShock = dir([DirectoryName '/' M_.fname '_inno*.mat']);

x=[];
yss = [];
logpo = [];
for j=1:length(fileParam),
    try,
    load([DirectoryName '/' M_.fname '_param' int2str(j) '.mat'],'stock','stock_logpo','stock_ys');
    x = [x; stock];
    yss = [yss; stock_ys];
    logpo = [logpo; stock_logpo];
    catch
    end
end

B = size(x,1);


inno = [];
ifilesmooth = 0;
ifileinno = 0;        
b = 0;
bsmooth=0;
binno = 0;
while b<B,
    b = b +1;
    if bsmooth<b,
        ifilesmooth = ifilesmooth + 1;
        load([DirectoryName '/' M_.fname '_smooth' int2str(ifilesmooth) '.mat'],'stock');
        smooth=stock;
        bsmooth0 = bsmooth;
        bsmooth = bsmooth + size(smooth,3);
    end
    
    if binno<b,
        ifileinno = ifileinno + 1;
        load([DirectoryName '/' M_.fname '_inno' int2str(ifileinno) '.mat'],'stock');
        inno=stock;
        binno0 = binno;
        binno = binno + size(inno,3);
    end
    set_all_parameters(x(b,:));
    [T,R,SteadyState]=dynare_resolve;
    if options_.loglinear
        alphahat=squeeze(smooth(oo_.dr.order_var,:,b-bsmooth0))-repmat(log(SteadyState(oo_.dr.order_var)),1,stock_gend);
    else
        alphahat=squeeze(smooth(oo_.dr.order_var,:,b-bsmooth0))-repmat(SteadyState(oo_.dr.order_var),1,stock_gend);
    end
    etahat=squeeze(inno(:,:,b-binno0));
    
    %           [alphahat,etahat,epsilonhat,alphatilde,SteadyState,trend_coeff,aK] = ...
    %           DsgeSmoother(x(b,:)',stock_gend,stock_data,[],0);
    
    aE=T*alphahat;
    ss_.T=T;
    ss_.R=R;
    
    ss_.a=alphahat(:,end);
    ss_.a1=alphahat(:,1);
    ss_.ss=SteadyState(oo_.dr.order_var);
    for j= (length(ss_.ss)+1):length(T),
        ss_.ss(j)=ss_.ss(find(T(j,:)));
    end
    ss_.etahat=etahat;
    ss_.alphahat=alphahat;
    % ss_.af=af;
    ss_.aE=aE;
    
    a1=ss_.a1;
    as=a1;
    a0 = a1-ss_.R*ss_.etahat(:,1);
    for j=1:size(ss_.etahat,2)-1,
        as(:,j+1)=ss_.T*as(:,j)+ss_.R*ss_.etahat(:,j+1);
    end
    gend=size(as,2);
    att=alphahat*0;
    inn=alphahat*0;
    deco=zeros(length(T), size(R,2), stock_gend);
    for j=2:stock_gend,
        att(:,j) = ss_.T^(j-1)*a1;
        inn(:,j) = ss_.R*ss_.etahat(:,j);
        if j>1,
            inn(:,j) = inn(:,j) +  ss_.T*inn(:,j-1);
        end
        for iexo=1:M_.exo_nbr,
            deco(:,iexo,j) = ss_.R(:,iexo)*ss_.etahat(iexo,j);
            if j>1,
                deco(:,iexo,j) = deco(:,iexo,j) +  ss_.T*deco(:,iexo,j-1);
            end
            
        end
    end
    if b==1,
        asM=as(oo_.dr.inv_order_var,:);
        attM=att(oo_.dr.inv_order_var,:);
        innM=inn(oo_.dr.inv_order_var,:);
        decoM=deco(oo_.dr.inv_order_var,:,:);
    else
        asM=asM+as(oo_.dr.inv_order_var,:);
        attM=attM+att(oo_.dr.inv_order_var,:);
        innM=innM+inn(oo_.dr.inv_order_var,:);
        decoM=decoM+deco(oo_.dr.inv_order_var,:,:);
    end
end

as = asM/B;
att = attM/B;
inn = innM/B;
deco = decoM/B;

for j=1:length(vname),
  clear sdec,
   indx = (strmatch(vname{j},M_.endo_names,'exact'));
  
  inno1 = inn(indx,kpoint1);
  inno2 = inn(indx,kpoint2);
  inno3 = inn(indx,kpoint3);
  inno4 = inn(indx,kpoint4);
  
  sdec1 = squeeze(deco(indx,:,kpoint1));
  sdec2 = squeeze(deco(indx,:,kpoint2));
  sdec3 = squeeze(deco(indx,:,kpoint3));
  sdec4 = squeeze(deco(indx,:,kpoint4));
  sdec0 = sdec1 + sdec2 + sdec3 + sdec4;
  
  for i=1:size(ex_names_,1),
    clear index,
    for ii=1:size(ex_names_,2),
      indbuf = strmatch(ex_names_{i,ii},M_.exo_names,'exact');
      if ~isempty(indbuf),
        index(ii) = indbuf;
      end
    end
%     sdec(i,:)=(sdec0(index,:));
    sdec(i,:)=sum(sdec0(index,:),1);
    sdec0(index,:)=0;
  end
  sdec=[sdec; sum(sdec0)];
  ipos=sdec>0;
  ineg=sdec<0;
  figure('name',vname{j})
  set(gcf,'position' ,[50 50 1100 850])
  hp=bar((sdec.*ipos)','stacked'); 
  for jjj=1:length(hp),
      set(get(hp(jjj),'children'),'Edgecolor',[0 0 0]);
  end
  %   set(get(hp(end),'children'),'Facecolor',[1 1 1]);
  hold on, hn=bar((sdec.*ineg)','stacked'); 
%   set(get(hn(end),'children'),'Facecolor',[1 1 1]);
  set(gca,'position',[0.05 0.3 0.9 0.65],'units','normalized')
  hleg=legend(leg,'interpreter','none','location','Best');
  set(hleg,'position',[0.35 0.05 0.3 0.2],'units','normalized')
  
  hleg2 = get(hleg,'children');
  cox=jet(length(leg)-1);
  cox=[cox; 1 1 1];
  colormap(cox);
%   set(get(hleg2(1),'children'),'Facecolor',[1 1 1]);
  hhh=get(gca,'children');
  for jjj=1:length(hhh), set(hhh(jjj),'Edgecolor',[0 0 0]), end,
  hold on, h1=plot(as(indx,kpoint1)+as(indx,kpoint2)+as(indx,kpoint3)+as(indx,kpoint4),'k:*');
  hold on, h1=plot(inno1+inno2+inno3+inno4,'k-d');
  set(h1,'MarkerFaceColor', 'k')  
  title(['year on year ',vname{j},' innovations'],'interpreter','none')
  set(gca,'Xtick',[1:4:length(kpoint1)+4])
  set(gca,'Xticklabel',TT(kpoint4(1):4:end))
  saveas(gcf,['post_shock_dec_yoy_',vname{j},'_2.fig'])
  eval(['print -dpdf post_shock_dec_yoy_',vname{j},'_2.pdf'])
  eval(['print -depsc2 post_shock_dec_yoy_',vname{j},'_2.eps'])
  set(gca,'xlim',[0 length(kpoint1)+4])
  
end
