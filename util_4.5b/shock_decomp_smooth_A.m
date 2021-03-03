function out = shock_decomp_smooth_A(xparam1, T0, ex_names_, leg, vname, no_initial_effect, texvname, nfrcst, init_names_, file_name)

global oo_ M_ options_
persistent flag_init

ngroups0 = size(ex_names_,1);
if size(ex_names_,2)>1 || ischar(ex_names_{1}),
    old_interface=1;
    ex_names_0 = cell([ngroups0 1]);
    for j=1:ngroups0,
        ex_names_0{j}=ex_names_(j,:);
    end
    ex_names_=ex_names_0;
elseif iscell(ex_names_{1})
    old_interface=0;
else
    error('Wrong formatting of of shocks/groups of shocks')
end    

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

if options_.TeX
    if flag_init
        fidTeX = fopen([M_.fname '_Shock_Decomposition_A.TeX'],'w+');
    else
        fidTeX = fopen([M_.fname '_Shock_Decomposition_A.TeX'],'a+');
    end
end

if nargin==1 && ~isempty(xparam1)
    set_all_parameters(xparam1);
end
ss_ = getSmootherInfo(M_,options_,oo_);
options_ = set_default_option(options_,'occbin_smoother',0);
if options_.occbin_smoother,
    %     load info1 regime_history
    regime_history = options_.occbin_regime_history;
    gend  = size(ss_.alphahat,2);
    nperiods_ = 30;
    irfshock =M_.exo_names;
    TT = repmat(ss_.T, [1 1 gend]);
    RR = repmat(ss_.R, [1 1 gend]);
    CC = zeros([size(ss_.R,1),gend ]);
    for tp = 1:gend-1,
        if any(regime_history(tp).regime1) || any(regime_history(tp).regime2),
            [~ ,~ ,~ , Tx, Rx, CONSTx, regime_history2] = mr_runsim_occbin_fn(0, options_.occbin_modnames, options_.occbin_constraints,irfshock,zeros(1,M_.exo_nbr),zeros(M_.endo_nbr,1),[],regime_history(tp), 1, nperiods_);
            TT(:,:,tp+1) = Tx(oo_.dr.order_var,oo_.dr.order_var);
            RR(:,:,tp+1) = Rx(oo_.dr.order_var,:);
            CC(:,tp+1) = CONSTx(oo_.dr.order_var);
        end
    end
end
if ls('data.mat')
    load data T;
end

TT=[0:0.25:ceil(size(ss_.etahat,2)/4)];
if exist('T')
    TT=T(options_.first_obs)+TT;
    %     [T(options_.first_obs:end):0.25:(max(T)+1)];
    % else
    %     TT=[1:size(ss_.etahat,2)+1];
end

if nargin<6 || isempty(no_initial_effect),
    no_initial_effect = 0;
end

conditional_decomp=0;
if nargin<2 || isempty(T0),
    t1decomp =1; % set it manually to change the initial point for plotting shock decomp
else
    if length(T0)>1,
        conditional_decomp=1;
        Tcond = T0(2);
        T0=T0(1);
        index_cond = find(TT==Tcond);
        ss_.a1=ss_.alphahat(:,index_cond);
        ss_.alphahat=ss_.alphahat(:,index_cond:end);
        ss_.etahat=ss_.etahat(:,index_cond:end);
        ss_.alphatt=ss_.alphatt(:,index_cond:end);
        ss_.aE=ss_.aE(:,index_cond:end);
        ss_.aEt=ss_.aEt(:,index_cond:end);
        TT=TT(index_cond:end);
    end
    t1decomp = max([1,find(TT==T0)]);
end

if (nargin<7) || (isempty(texvname) )
    for j=1:length(vname),
        iendo(j) =  strmatch(vname{j},M_.endo_names,'exact');
        texvname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end
if (nargin<8) || (isempty(nfrcst) )
    nfrcst=0;
end
if (nargin<9) || isempty(init_names_),
    init_names_=[];
    indx_init=[];
end
for j=1:length(init_names_)
    indx_init(j)=strmatch(init_names_{j},M_.endo_names(oo_.dr.order_var,:),'exact');
end
% re-compute smoothed variables from T, R and  smoothed shocks etahat
a1=ss_.a1;
as=a1;
as0=zeros(length(a1),1);
as0(indx_init)=a1(indx_init);
a0 = a1-ss_.R*ss_.etahat(:,1);
ss_.etahat=[ss_.etahat zeros(size(ss_.etahat,1),nfrcst)];
for j=1:size(ss_.etahat,2)-1,
    if options_.occbin_smoother,
        TM = TT(:,:,j+1);
        RM = RR(:,:,j+1);
        CONST = CC(:,j+1);
    else
        TM = ss_.T;
        RM= ss_.R;
        CONST = a1*0;
    end
    as(:,j+1)=TM*as(:,j)+RM*ss_.etahat(:,j+1)+CONST;
    as0(:,j+1)=TM*as0(:,j)+CONST;
end
gend=size(as,2);
att=zeros(M_.endo_nbr,gend);
inn=zeros(M_.endo_nbr,gend);
deco=zeros(M_.endo_nbr,M_.exo_nbr,gend);
att(:,1)=a1;
for j=2:gend,
    if options_.occbin_smoother,
        TM = TT(:,:,j);
        RM = RR(:,:,j);
        CONST = CC(:,j);
    else
        TM = ss_.T;
        RM= ss_.R;
        CONST = a1*0;
    end
    att(:,j) = TM*att(:,j-1)+CONST;
    inn(:,j) = RM*ss_.etahat(:,j);
    if j>1,
        inn(:,j) = inn(:,j) +  TM*inn(:,j-1); %+CONST;
    end
    for iexo=1:M_.exo_nbr,
        deco(:,iexo,j) = RM(:,iexo)*ss_.etahat(iexo,j);
        if j>1,
            deco(:,iexo,j) = deco(:,iexo,j) +  TM*deco(:,iexo,j-1);%+CONST;
        end
        
    end
end
as=as(oo_.dr.inv_order_var,:);
as0=as0(oo_.dr.inv_order_var,:);
att=att(oo_.dr.inv_order_var,:);
inn=inn(oo_.dr.inv_order_var,:);
deco=deco(oo_.dr.inv_order_var,:,:);


% kpoint = t1decomp:options_.nobs;
kpoint = t1decomp:gend;
% kpoint3 = kpoint4-1;
% kpoint2 = kpoint4-2;
% kpoint1 = kpoint4-3;
% kpoint1 = 62:101;

ngroups = ngroups0+1+no_initial_effect;
ncol=3;
nrow=ceil(ngroups/ncol);
% nrow=round(sqrt(ngroups));
% if nrow^2 < ngroups,
%     ncol=nrow+1;
% else
%     ncol=nrow;
% end

for j=1:length(vname),
    clear sdec  sdec2,
    
    indx = (strmatch(vname{j},M_.endo_names,'exact'));
    
    inno = inn(indx,kpoint);
    %   inno2 = inn(indx,kpoint2);
    %   inno3 = inn(indx,kpoint3);
    %   inno4 = inn(indx,kpoint4);
    
    sdec0 = squeeze(deco(indx,:,kpoint));
    %   sdec2 = squeeze(deco(indx,:,kpoint2));
    %   sdec3 = squeeze(deco(indx,:,kpoint3));
    %   sdec4 = squeeze(deco(indx,:,kpoint4));
    
    for i=1:ngroups0,
        clear index,
        for ii=1:size(ex_names_{i},2),
            indbuf = strmatch(ex_names_{i}{ii},M_.exo_names,'exact');
            if ~isempty(indbuf),
                index(ii) = indbuf;
            elseif ~isempty(ex_names_{i}{ii}),
                error(['Shock name ',ex_names_{i}{ii}, ' not found.' ]);
            end
        end
        %     sdec(i,:)=(sdec0(index,:));
        sdec(i,:)=sum(sdec0(index,:),1);
        sdec0(index,:)=0;
    end
    if ~isempty(indx_init)
        sdec=[sdec; as0(indx,kpoint)];
    end
    sdec2=[sdec; as(indx,kpoint)-inno+sum(sdec0)-as0(indx,kpoint)];
    sdec=[sdec; sum(sdec0)];
    nbcmpts=size(sdec,1);
    
    
    if flagmap==0
        func = @(x) colorspace('RGB->Lab',x);
        MAP = distinguishable_colors(nbcmpts,'w',func);
        MAP(end,:) = [0.7 0.7 0.7];
    else
        colours=leg(:,2);
        
        for i=1:length(colours); MAP(i,:)=rgb(colours{i,:}); end
    end
    
    
    leg0=leg(:,1);
    
    ipos=sdec>0;
    ineg=sdec<0;
    ipos2=sdec2>0;
    ineg2=sdec2<0;
    if no_initial_effect,
        
        h = dyn_figure(options_.nodisplay,'name',vname{j}, 'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','landscape','renderermode','auto');
        colormap(MAP);
        set(gcf,'position' ,[50 50 1100 850])
        hp=bar((sdec(:,4:4:end).*ipos(:,4:4:end))','stacked','EdgeColor',[0 0 0]);
        shading faceted;
        %   set(get(hp(end),'children'),'Facecolor',[1 1 1]);
        hold on, hn=bar((sdec(:,4:4:end).*ineg(:,4:4:end))','stacked');
        shading faceted;
        %   set(get(hn(end),'children'),'Facecolor',[1 1 1]);
        
        set(gca,'position',[0.05 0.23 0.95 0.72],'units','normalized')
        hleg=legend(leg0,'interpreter','none','location','BestOutside');
        shading faceted;
%        set(hleg,'position',[0.35 0.05 0.3 0.3],'units','normalized')
        
        hleg2 = get(hleg,'children');
        child = get(hleg2(1:2:end),'children');
        for ichild = 1:length(child),
            set(child{ichild},'Edgecolor',[0 0 0]);
        end
        if ~no_initial_effect,
            hold on, h1=plot(as(indx,kpoint(:,4:4:end)),'k-d');
            set(h1,'MarkerFaceColor', 'k')
        end
        hold on, h1=plot(inno(4:4:end),'k:*');
        %   if nargin==7
        %       title(['quarter on quarter ',texvname{j},' innovations'],'interpreter','tex')
        %   else
        tit1 = ['Annual ',vname{j},' shock decomposition'];
        title(tit1,'interpreter','none')
        %   end
        SS0=[1:4:length(kpoint)];
        SS1=1:length(SS0);
        set(gca,'Xtick',SS1);
        set(gca,'Xticklabel',floor(TT(SS0+kpoint(1))));
        set(gca,'xlim',[0 length(kpoint(4:4:end))+1])
        a=axis;
        if a(3)+get_mean(vname{j})<0 && a(4)+get_mean(vname{j})>0,
            plot(a(1:2),get_mean(vname{j})*[-1 -1],'k--','linewidth',1)
            ytick=get(gca,'ytick');
            ytick1=ytick-get_mean(vname{j});
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
        set(gca,'yticklabel',num2str(get(gca,'ytick')'+get_mean(vname{j}),'%4.2g'))
        dyn_saveas(gcf,[M_.fname,'_shock_dec_AoA_',vname{j},'_vs_',num2str(ngroups0),'Shocks',file_name],options_.nodisplay,options_.graph_format);
        
    else
        hinit = dyn_figure(options_.nodisplay,'name',vname{j}, 'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','landscape','renderermode','auto');
        colormap(MAP);
        set(gcf,'position' ,[50 50 1100 850])
        hp=bar((sdec2(:,4:4:end).*ipos2(:,4:4:end))','stacked','EdgeColor',[0 0 0]);
        shading faceted;
        %   set(get(hp(end),'children'),'Facecolor',[1 1 1]);
        hold on, hn=bar((sdec2(:,4:4:end).*ineg2(:,4:4:end))','stacked');
        shading faceted;
        %   set(get(hn(end),'children'),'Facecolor',[1 1 1]);
        
        set(gca,'position',[0.05 0.23 0.95 0.72],'units','normalized')
        hleg=legend(leg0,'interpreter','none','location','BestOutside');
        shading faceted;
%         set(hleg,'position',[0.8 0.38 0.19 0.57],'units','normalized')
        
        hleg2 = get(hleg,'children');
        child = get(hleg2(1:2:end),'children');
        for ichild = 1:length(child),
            set(child{ichild},'Edgecolor',[0 0 0]);
        end
        hold on, h1=plot(as(indx,kpoint(:,4:4:end)),'k-d');
        set(h1,'MarkerFaceColor', 'k')
        %   if nargin==7
        %       title(['quarter on quarter ',texvname{j}],'interpreter','tex')
        %   else
        tit2 = ['Annual ',vname{j}];
        title(tit2,'interpreter','none')
        %   end
        SS0=[1:4:length(kpoint)];
        SS1=1:length(SS0);
        set(gca,'Xtick',SS1);
        set(gca,'Xticklabel',floor(TT(SS0+kpoint(1))));
        set(gca,'xlim',[0 length(kpoint(4:4:end))+1])
        a=axis;
        if a(3)+get_mean(vname{j})<0 && a(4)+get_mean(vname{j})>0,
            plot(a(1:2),get_mean(vname{j})*[-1 -1],'k--','linewidth',1)
            ytick=get(gca,'ytick');
            ytick1=ytick-get_mean(vname{j});
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
        set(gca,'yticklabel',num2str(get(gca,'ytick')'+get_mean(vname{j}),'%4.2g'))
        dyn_saveas(gcf,[M_.fname,'_shock_dec_AoA_',vname{j},'_vs_',num2str(ngroups0),'Shocks','_Init',file_name],options_.nodisplay,options_.graph_format)
    end
    hdetail = dyn_figure(options_.nodisplay,'name',vname{j}, 'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','portrait','renderermode','auto','position',[200 100 650 850]);
    a0=a*0;
    a0(3)=inf;
    a0(4)=-inf;
    for i=1:ngroups0+1+no_initial_effect,
        if i==ngroups0+2;
            distr=sum(sdec2)-sum(sdec);
            distr(2,:)=sum(sdec2)-distr;
        else
            if no_initial_effect
                distr=sdec(i,:);
            else
                distr=sdec2(i,:);
            end
            distr(2,:)=sum(sdec2)-distr;
        end
        ipos=distr>0;
        ineg=distr<0;
        subplot(nrow,ncol,i), set(gca,'box','on')
        hbar = bar((distr(:,4:4:end).*ipos(:,4:4:end))','stacked');
        colormap([0.15 0.15 0.15;0.85 0.85 0.85]),
        %         shading faceted;
        set(hbar,'edgecolor','flat');
        hold on,
        hbar = bar((distr(:,4:4:end).*ineg(:,4:4:end))','stacked');
        colormap([0.15 0.15 0.15;0.85 0.85 0.85]),
        %         shading faceted;
        set(hbar,'edgecolor','flat');
        %         bar(distr','facecolor',MAP(end,:),'edgecolor',MAP(end,:)),
        if i==ngroups0+2;
            title('Initial condition')
        else
            title(leg0{i}),
        end
        axis tight;
        a=axis;
        if length(get(gca,'Xtick'))>5
            SS0=[1:16:length(kpoint)];
            SS1=1:4:length(kpoint)/4;
            set(gca,'Xtick',SS1);
            set(gca,'Xticklabel',floor(TT(SS0+kpoint(1))));
        else
            SS0=[1:4:length(kpoint)];
            set(gca,'Xticklabel',floor(TT(SS0(get(gca,'Xtick'))+kpoint(1))));
        end
        set(gca,'xlim',[0 length(kpoint(4:4:end))+1])
        a0(3)=min(a(3),a0(3));
        a0(4)=max(a(4),a0(4));
        set(gca,'ylim',a0(3:4))
        hold on, h1=plot(as(indx,kpoint(4:4:end)),'k-');
        set(h1,'MarkerFaceColor', 'k','linewidth',2)
    end
    for i=1:ngroups0+1+no_initial_effect,
        subplot(nrow,ncol,i),
        set(gca,'ylim',a0(3:4))
    end
    %     nodisplay0= options_.nodisplay;
    %     options_.nodisplay = 0;
    %     if nodisplay0, set(gcf,'visible','on'), end
    if ~no_initial_effect,
        dyn_saveas(gcf,[M_.fname,'_shock_dec_AoA_',vname{j},'_vs_',num2str(ngroups0),'Shocks','_Init_Detail',file_name],options_.nodisplay,options_.graph_format)
    else
        %     if nodisplay0, set(gcf,'visible','off'), end
        %     options_.nodisplay = nodisplay0;
        %     go=subplot(nrow,ncol,i);
        %     delete(go);
        %     go=subplot(nrow,ncol,i-1);
        %     delete(go);
        %     subplot(nrow,ncol,i-1), set(gca,'box','on')
        %     distr=sdec2(end,:);
        %     distr(2,:)=sum(sdec2)-distr;
        %     hbar = bar((distr(:,4:4:end).*(distr(:,4:4:end)>0))','stacked');
        %     colormap([0.15 0.15 0.15;0.85 0.85 0.85]),
        %     %         shading faceted;
        %     set(hbar,'edgecolor','flat');
        %     hold on,
        %     hbar = bar((distr(:,4:4:end).*(distr(:,4:4:end)<0))','stacked');
        %     colormap([0.15 0.15 0.15;0.85 0.85 0.85]),
        %     %         shading faceted;
        %     set(hbar,'edgecolor','flat');
        %     %         bar(distr','facecolor',MAP(end,:),'edgecolor',MAP(end,:)),
        %     title(leg0{i-1}),
        %     axis tight;
        %     a=axis;
        %     SS0=[1:16:length(kpoint)];
        %     SS1=1:4:length(kpoint)/4;
        %     set(gca,'Xtick',SS1);
        %     set(gca,'Xticklabel',floor(TT(SS0+kpoint(1))));
        %     set(gca,'xlim',[0 length(kpoint(4:4:end))+1])
        %     a0(3)=min(a(3),a0(3));
        %     a0(4)=max(a(4),a0(4));
        %     set(gca,'ylim',a0(3:4))
        %     hold on, h1=plot(as(indx,kpoint(4:4:end)),'k-');
        %     set(h1,'MarkerFaceColor', 'k','linewidth',2)
        dyn_saveas(gcf,[M_.fname,'_shock_dec_AoA_',vname{j},'_vs_',num2str(ngroups0),'Shocks','_Detail',file_name],options_.nodisplay,options_.graph_format)
    end
    tempshockdecom=['AoA_',vname{j},'_vs_',num2str(ngroups0),'_Shocks',file_name];
    eval([tempshockdecom,'.Decomp','=sdec(:,4:4:end);']);
    eval([tempshockdecom,'.Init_Decomp','=sdec2(:,4:4:end);']);
    eval([tempshockdecom,'.Smooth','=inno(4:4:end);']);
    eval([tempshockdecom,'.Init_Smooth','=as(indx,kpoint(4:4:end));']);
    eval([tempshockdecom,'.Shocks','=ex_names_;']);
    eval([tempshockdecom,'.Names','=leg;']);
    eval([tempshockdecom,'.InitialTime','=T0;']);
    TagName = [M_.fname,'_AoA_ShockDecomp_vs_',num2str(ngroups0),'Shocks',file_name];
    eval([tempshockdecom,'.TagName=TagName;']);
    SheetName = [vname{j}];
    eval([tempshockdecom,'.SheetName=SheetName;']);
    
    
    if nargout>0,
        out.shocks=ex_names_;
        out.names=leg;
        eval(['out.decomp',vname{j},'=sdec(:,4:4:end);']);
        eval(['out.decomp2',vname{j},'=sdec2(:,4:4:end);']);
        eval(['out.inno',vname{j},'=inno1(:,4:4:end));']);
        eval(['out.smooth',vname{j},'=as(indx,kpoint1(4:4:end));']);
        
    end
    
    if options_.TeX
        
        if no_initial_effect,
            tempfilename=[M_.fname,'_shock_dec_AoA_',vname{j},'_vs_',num2str(ngroups0),'Shocks',file_name];
            
            [a,b,c]=fileparts(tempfilename);
            
            
            % % %             [a,b,c]=fileparts(shockfilename);
            % % %             %     quest3hlmr_shock_dec_AoAt_E_LCY_Init
            % % %             pos=findstr(b,'_');
            % % %             varname=b(pos(end-2)+1:pos(end)-1);
            % % %
            % % %             varpos= strmatch(varname, M_.endo_names, 'exact');
            % % %             texname=M_.endo_names_tex(varpos,:);
            % % %
            
            % TeX eps loader file
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.8][0]{%s}\n',tit1,['Annual $' deblank(texvname{j}) '$ innovations']);
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,['\\includegraphics[width=0.65\\textwidth,angle=-90] {' b '} \n']);
            fprintf(fidTeX,'\\caption{Shock Decomposition Annual $%s$ innovations}',texvname{j});
            fprintf(fidTeX,'\\label{Fig:Shock Decomposition:%s}\n',int2str(j));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
            
            tempfilename=[M_.fname,'_shock_dec_AoA_',vname{j},'_vs_',num2str(ngroups0),'Shocks','_Detail',file_name];
            
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
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.8][0]{%s}\n',tit2,['Annual $' deblank(texvname{j}) '$']);
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,['\\includegraphics[width=0.65\\textwidth,angle=-90] {' b '} \n']);
            fprintf(fidTeX,'\\caption{Shock Decomposition Annual $%s$}',texvname{j});
            fprintf(fidTeX,'\\label{Fig:Shock Decomposition Detail:%s}\n',int2str(j));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
            
        else
            
            tempfilename=[M_.fname,'_shock_dec_AoA_',vname{j},'_vs_',num2str(ngroups0),'Shocks','_Init',file_name];
            
            [a,b,c]=fileparts(tempfilename);
            
            
            % % %             [a,b,c]=fileparts(shockfilename);
            % % %             %     quest3hlmr_shock_dec_AoAt_E_LCY_Init
            % % %             pos=findstr(b,'_');
            % % %             varname=b(pos(end-2)+1:pos(end)-1);
            % % %
            % % %             varpos= strmatch(varname, M_.endo_names, 'exact');
            % % %             texname=M_.endo_names_tex(varpos,:);
            % % %
            
            % TeX eps loader file
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.8][0]{%s}\n',tit2,['Annual $' deblank(texvname{j}) '$']);
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,['\\includegraphics[width=0.65\\textwidth,angle=-90] {' b '} \n']);
            fprintf(fidTeX,'\\caption{Shock Decomposition Annual $%s$}',texvname{j});
            fprintf(fidTeX,'\\label{Fig:Shock Decomposition:%s}\n',int2str(j));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
            
            
            tempfilename=[M_.fname,'_shock_dec_AoA_',vname{j},'_vs_',num2str(ngroups0),'Shocks','_Init_Detail',file_name];
            
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
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.8][0]{%s}\n',tit2,['Annual $' deblank(texvname{j}) '$']);
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,['\\includegraphics[width=0.65\\textwidth,angle=-90] {' b '} \n']);
            fprintf(fidTeX,'\\caption{Shock Decomposition annual $%s$}',texvname{j});
            fprintf(fidTeX,'\\label{Fig:Shock Decomposition Init Detail:%s}\n',int2str(j));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
        end
    end
    if (flag_init && j==1) || isempty(ls([M_.fname,'_AoA_ShockDecomp.mat']))
        %         delete  *shock_dec*AoA*.fig
        %         delete  *shock_dec*AoA*.eps
        save([M_.fname '_AoA_ShockDecomp'],tempshockdecom);
    else
        save([M_.fname '_AoA_ShockDecomp'],tempshockdecom,'-append');
    end
end


if options_.TeX
    fclose(fidTeX);
end