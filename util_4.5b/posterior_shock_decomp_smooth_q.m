function posterior_shock_decomp_smooth_q(M_, options_, oo_, estim_params_, T0, ex_names_, leg, vname, texvname, file_name)

persistent flag_init

if nargin<10 || isempty(file_name)
    file_name='';
end

if isempty(flag_init)
    flag_init=1;
else
    flag_init=0;
end

if size(leg,2)==1,
    flagmap=0;
else
    flagmap=1;
end

DirectoryName = CheckPath('metropolis',M_.dname);
OutDirectoryName = CheckPath('Output',M_.dname);

if options_.TeX
    if flag_init
        fidTeX = fopen([OutDirectoryName,filesep,M_.fname '_Posterior_Shock_Decomposition_q.TeX'],'w+');
    else
        fidTeX = fopen([OutDirectoryName,filesep,M_.fname '_Posterior_Shock_Decomposition_q.TeX'],'a+');
    end
end

load data T;

if (nargin<9) || (isempty(texvname) )
    for j=1:length(vname),
        iendo(j) =  strmatch(vname{j},M_.endo_names,'exact');
        texvname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end


for j=1:length(vname),  
    indx(j) = find(strcmp(vname{j},cellstr(M_.endo_names)));
end
for i=1:size(ex_names_,1),
    for ii=1:size(ex_names_,2),
        indbuf = strmatch(ex_names_{i,ii},M_.exo_names,'exact');
        if ~isempty(indbuf),
            index{i}(ii) = indbuf;
        end
    end
end

smooth_file_list = dir([DirectoryName,filesep,'*_smooth*.mat']);
inno_file_list = dir([DirectoryName,filesep,'*_inno*.mat']);
tmp_file_list = dir([DirectoryName,filesep,'*_param*.mat']);
jfile=0;
for j=1:length(tmp_file_list),
    if isempty(strfind(tmp_file_list(j).name,'irf')),
        jfile=jfile+1;
        param_file_list(jfile)=tmp_file_list(j);
    end
end
clear tmp_file_list jfile,
jsmooth_file=1;
jinno_file=1;
jsmooth=0;
jinno=0;
stock_inno = load([DirectoryName,filesep,inno_file_list(1).name],'stock');
load([DirectoryName,filesep,smooth_file_list(1).name],'stock')
gend=size(stock,2);

TT=[0:0.25:ceil(gend/4)];
if exist('T')
    TT=T(options_.first_obs)+TT;
    %     [T(options_.first_obs:end):0.25:(max(T)+1)];
    % else
    %     TT=[1:size(ss_.etahat,2)+1];
end

if nargin<5 || isempty(T0),
    t1decomp =25; % set it manually to change the initial point for plotting shock decomp
else
    t1decomp = max(1,find(TT==T0));
end

kpoint = t1decomp:gend;

B=0;
SS=[];
for j=1:length(param_file_list),
    load([DirectoryName,filesep,param_file_list(j).name],'stock_ys')
    stock_params=load([DirectoryName,filesep,param_file_list(j).name],'stock');
    stock_decomp=zeros(length(vname),size(ex_names_,1)+1,gend,size(stock_params.stock,1));
    B=B+size(stock_params.stock,1);
    SS=[SS;stock_ys(1:size(stock_params.stock,1),indx)];
    for i=1:size(stock_params.stock,1),
        jsmooth = jsmooth+1;
        if jsmooth>size(stock,3),
            jsmooth=1;
%             stock=stock_decomp;
%             save([DirectoryName,filesep,M_.fname,'_decomp',int2str(jsmooth_file)],'stock')
            jsmooth_file = jsmooth_file+1;
            load([DirectoryName,filesep,M_.fname,'_smooth',int2str(jsmooth_file),'.mat'],'stock')
%             stock_decomp=zeros(length(vname),size(ex_names_,1)+1,gend,size(stock,3));
        end
        ahat = squeeze(stock(:,:,jsmooth));
        jinno = jinno+1;
        if jinno>size(stock_inno.stock,3),
            jinno=1;
            jinno_file = jinno_file+1;
            stock_inno = load([DirectoryName,filesep,M_.fname,'_inno',int2str(jinno_file),'.mat'],'stock');
        end
        ehat = squeeze(stock_inno.stock(:,:,jinno));
        M_ = set_all_parameters(stock_params.stock(i,:)',estim_params_,M_);
        [T,R,SteadyState,info] = dynare_resolve(M_, options_, oo_);
        if  info,
            options_.aim_solver=1;
            [T,R,SteadyState,info] = dynare_resolve(M_, options_, oo_);
            options_.aim_solver=0;
        end
        as=ahat-repmat(SteadyState,1,gend);
        ahat=ahat(oo_.dr.order_var,:)-repmat(SteadyState(oo_.dr.order_var),1,gend);
%         ahat=ahat(oo_.dr.order_var,:)-repmat(stock_ys(i,oo_.dr.order_var)',1,size(ahat,2));

        a1=ahat(:,1);
%         as=a1;
%         for jx=1:gend-1,
%             as(:,jx+1)=T*as(:,jx)+R*ehat(:,jx+1);
%         end
        inn=zeros(M_.endo_nbr,gend);
        att=zeros(M_.endo_nbr,gend);
        deco=zeros(M_.endo_nbr,M_.exo_nbr,gend);
        att(:,1)=a1;
        for jx=2:gend,
            att(:,jx) = T*att(:,jx-1);
            inn(:,jx) = R*ehat(:,jx);
            if jx>2,
                inn(:,jx) = inn(:,jx) +  T*inn(:,jx-1);
            end
            for iexo=1:M_.exo_nbr,
                deco(:,iexo,jx) = R(:,iexo)*ehat(iexo,jx);
                if jx>1,
                    deco(:,iexo,jx) = deco(:,iexo,jx) +  T*deco(:,iexo,jx-1);
                end
                
            end
        end
%         as=as(oo_.dr.inv_order_var,:);
        inn=inn(oo_.dr.inv_order_var,:);
        deco=deco(oo_.dr.inv_order_var,:,:);
        inno = inn(indx,kpoint);
        
        sdec0 = deco(indx,:,kpoint);
        
        for jx=1:length(vname),
            for ix=1:size(ex_names_,1),
                sdec(jx,ix,:)=sum(sdec0(jx,index{ix},:),2);
                sdec0(jx,index{ix},:)=0;
            end
            sdec2(jx,:,:)=[squeeze(sdec(jx,:,:)); as(indx(jx),kpoint)-inno(jx,:)+sum(squeeze(sdec0(jx,:,:)),1)];
        end
        stock_decomp(:,:,:,i) = sdec2;
    end
    stock0=stock;
    stock=stock_decomp;
    save([DirectoryName,filesep,M_.fname,'_decomp_vs_',num2str(size(ex_names_,1)),'Shocks',file_name,'_',int2str(j)],'stock')
    stock=stock0;
    clear stock0;
end
stock1 = zeros(length(vname),size(ex_names_,1)+1,gend,B);
k = 0;
for file = 1:length(param_file_list)
    load([DirectoryName,filesep,M_.fname,'_decomp_vs_',num2str(size(ex_names_,1)),'Shocks',file_name,'_',int2str(file)],'stock')
    k = k(end)+(1:size(stock,4));
    stock1(:,:,:,k) = stock;
end
clear stock
SSMean = mean(SS);
for j = 1:length(vname),
    stock2=sum(stock1,2);
    for is = 1:gend
        [YMean(j,is),YMedian(j,is),YVar(j,is),YHPD(:,j,is),YDistrib(:,j,is)] = ...
            posterior_moments(squeeze(stock2(j,1,is,:)),0,options_.mh_conf_sig);
    end
for i=1:size(ex_names_,1)+1,
    for is = 1:gend
        [Mean(j,i,is),Median(j,i,is),Var(j,i,is),HPD(:,j,i,is),Distrib(:,j,i,is)] = ...
            posterior_moments(squeeze(stock1(j,i,is,:)),0,options_.mh_conf_sig);
    end
end
end

ngroups = size(ex_names_,1)+1;
nrow=round(sqrt(ngroups));
if nrow^2 < ngroups,
    ncol=nrow+1;
else
    ncol=nrow;
end

for j=1:length(vname),
    nbcmpts=size(Mean,2);
       
    if flagmap==0
        MAP=CreateColorMap(nbcmpts);
        
        MAP(end,:) = [0.9 0.9 0.9]; 
    else
        colours=leg(:,2);
        
        for i=1:length(colours); MAP(i,:)=rgb(colours{i,:}); end
    end
    
    
    leg0=leg(:,1);
    
    sdec2=squeeze(Mean(j,:,:));
    as=sum(sdec2,1);
    ipos2=sdec2>0;
    ineg2=sdec2<0;
        
    hinit = dyn_figure(options_.nodisplay,'name',vname{j}, 'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','landscape','renderermode','auto');
    colormap(MAP);
    set(gcf,'position' ,[50 50 1100 850])
    hp=bar((sdec2.*ipos2)','stacked','EdgeColor',[0 0 0]);
    shading faceted;
    %   set(get(hp(end),'children'),'Facecolor',[1 1 1]);
    hold on, hn=bar((sdec2.*ineg2)','stacked');
    shading faceted;
    %   set(get(hn(end),'children'),'Facecolor',[1 1 1]);
    
    set(gca,'position',[0.05 0.4 0.9 0.55],'units','normalized')
    hleg=legend(leg0,'interpreter','none','location','Best');
    shading faceted;
    set(hleg,'position',[0.35 0.05 0.3 0.3],'units','normalized')
    
    hleg2 = get(hleg,'children');
    child = get(hleg2(1:2:end),'children');
    for ichild = 1:length(child),
        set(child{ichild},'Edgecolor',[0 0 0]);
    end
    hold on,
    plot(squeeze(YHPD(:,j,:))',':k')
    h1=plot(as(kpoint),'k-d');
    set(h1,'MarkerFaceColor', 'k')
    %   if nargin==7
    %       title(['quarter on quarter ',texvname{j}],'interpreter','tex')
    %   else
    tit2 = ['quarter to quarter ',vname{j}];
    title(tit2,'interpreter','none')
    %   end
    set(gca,'Xtick',[1:4:length(kpoint)])
    set(gca,'Xticklabel',TT(kpoint(1):4:end))
    set(gca,'xlim',[0 length(kpoint)+4])
    a=axis;
    if a(3)+SSMean(j)<0 && a(4)+SSMean(j)>0,
        plot(a(1:2),SSMean(j)*[-1 -1],'k--','linewidth',1)
        ytick=get(gca,'ytick');
        ytick1=ytick-SSMean(j);
        ind1=min(find(ytick1>=a(3)));
        ind2=max(find(ytick1<=a(4)));
        dytick=ytick(2)-ytick(1);      
        if ind1>1,
            ytick1  = [ytick1(ind1:end) ytick1(end)+dytick:dytick:a(4)];
        elseif ind2<length(ytick)
            ytick1= [sort(ytick1(1)-dytick:-dytick:a(3)) ytick1(1:ind2)];
        end
        set(gca,'ytick',ytick1),
    end
    set(gca,'yticklabel',num2str(get(gca,'ytick')'+SSMean(j),'%4.2g'))
    dyn_saveas(gcf,[OutDirectoryName,filesep,M_.fname,'_posterior_shock_dec_qoq_',vname{j},'_vs_',num2str(size(ex_names_,1)),'Shocks',file_name,'_Init'],options_.nodisplay,options_.graph_format)

    hinit = dyn_figure(options_.nodisplay,'name',vname{j}, 'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','landscape','renderermode','auto');
    for i=1:size(ex_names_,1)+1,
        distr=squeeze(Distrib(:,j,i,:));
        subplot(nrow,ncol,i), plot(distr','color',MAP(i,:)), title(leg0{i}),
        set(gca,'Xtick',[1:16:length(kpoint)])
        set(gca,'Xticklabel',TT(kpoint(1):16:end))
        set(gca,'xlim',[0 length(kpoint)+4])
        hold all, plot(sdec2(i,:),'k')
    end
    dyn_saveas(gcf,[OutDirectoryName,filesep,M_.fname,'_posterior_shock_dec_qoq_',vname{j},'_vs_',num2str(size(ex_names_,1)),'Shocks',file_name,'_Init_Detail'],options_.nodisplay,options_.graph_format)
    if options_.TeX
        
        tempfilename=[OutDirectoryName,filesep,M_.fname,'_posterior_shock_dec_qoq_',vname{j},'_vs_',num2str(size(ex_names_,1)),'Shocks',file_name,'_Init'];
        
        [a,b,c]=fileparts(tempfilename);
        
        
        % % %             [a,b,c]=fileparts(shockfilename);
        % % %             %     quest3hlmr_shock_dec_qoqt_E_LCY_Init
        % % %             pos=findstr(b,'_');
        % % %             varname=b(pos(end-2)+1:pos(end)-1);
        % % %
        % % %             varpos= strmatch(varname, M_.endo_names, 'exact');
        % % %             texname=M_.endo_names_tex(varpos,:);
        % % %
        
        % TeX eps loader file
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.8][0]{%s}\n',tit2,['quarter on quarter $' deblank(texvname{j}) '$']);
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,['\\includegraphics[width=0.65\\textwidth,angle=-90] {' b '} \n']);
        fprintf(fidTeX,'\\caption{Shock Decomposition quarter on quarter $%s$}',texvname{j});
        fprintf(fidTeX,'\\label{Fig:Shock Decomposition:%s}\n',int2str(j));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
        
    end
    
end

if options_.TeX
    fclose(fidTeX);
end
