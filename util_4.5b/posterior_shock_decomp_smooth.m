function posterior_shock_decomp_smooth(M_,oo_,options_,varlist)
% function posterior_shock_decomp_smooth(M_,oo_,options_,varlist)
% XXX [TO COMPLETE]
%
% INPUTS
%    M_:          [structure]  Definition of the model
%    oo_:         [structure]  Storage of results
%    options_:    [structure]  Options
%    varlist:     [char]       List of variables
%
% SPECIAL REQUIREMENTS
%    none
%
% Copyright (C) 2016-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.


%% Declarations and set-up
% Tests whether we already entered the function (NO=>flag_init=1; 0 otherwise)
persistent flag_init
global estim_params_;

if isempty(flag_init)
    flag_init=1;
else
    flag_init=0;
end

% File name
if(options_.use_shock_groups)
    file_name = sprintf('_%s',options_.use_shock_groups);
else
    file_name='';
end

% The group colors are not given (=0)
flagmap=0;

% Reads the options_
type    =   options_.shock_decomp.type;
%T0      =   options_.shock_decomp.init_state; % Warning, also we might consider options_.initial_date
T0      =   options_.initial_date.double();

% Set-up the figure file names if provided
fig_names='';
if ~isempty(options_.plot_shock_decomp.fig_name)
	fig_names=['_' options_.plot_shock_decomp.fig_name];
end


%
DirectoryName = CheckPath('metropolis',M_.dname);
OutDirectoryName = CheckPath('graphs',M_.dname);


% Sets and opens the TeX document
switch type
    case 'qoq' 
        if options_.TeX
            if flag_init
                fidTeX = fopen([OutDirectoryName,filesep,M_.fname '_Posterior_Shock_Decomposition_q.TeX'],'w+');
            else
                fidTeX = fopen([OutDirectoryName,filesep,M_.fname '_Posterior_Shock_Decomposition_q.TeX'],'a+');
            end
        end
    case 'yoy'
        if options_.TeX
            if flag_init
                fidTeX = fopen([OutDirectoryName,filesep,M_.fname '_Posterior_Shock_Decomposition_yoy.TeX'],'w+');
            else
                fidTeX = fopen([OutDirectoryName,filesep,M_.fname '_Posterior_Shock_Decomposition_yoy.TeX'],'a+');
            end
        end                          
    case ''
        if options_.TeX
            if flag_init
                fidTeX = fopen([OutDirectoryName,filesep,M_.fname '_Posterior_Shock_Decomposition.TeX'],'w+');
            else
                fidTeX = fopen([OutDirectoryName,filesep,M_.fname '_Posterior_Shock_Decomposition.TeX'],'a+');
            end
        end                          

    otherwise
        error('posterior_shock_decomp_smooth:: Wrong type (options_.shock_decomp.type). Should be "qoq", "yoy" or "aoa"')
end   
    

%%
% number of variables
endo_nbr = M_.endo_nbr;
% number of shocks
nshocks = M_.exo_nbr;
% Retrieves the shock names as listed in options_
% Remark that 'old_interface' is not used.
z = oo_.shock_decomposition;
gend = size(z,3);
if options_.use_shock_groups
    shock_groups = M_.shock_groups.(options_.use_shock_groups);
    shock_ind = fieldnames(shock_groups);
    ngroups = length(shock_ind);
    fig_names=[fig_names ' group ' options_.use_shock_groups];
    shock_names = shock_ind;
    for i=1:ngroups,
       shock_names{i} = (shock_groups.(shock_ind{i}).label);
    end
    zz = zeros(endo_nbr,ngroups+2,gend);
    kcum=[];
    for i=1:ngroups
        for j = shock_groups.(shock_ind{i}).shocks
            k = find(strcmp(j,cellstr(M_.exo_names)));
            zz(:,i,:) = zz(:,i,:) + z(:,k,:);
            z(:,k,:) = 0;
            kcum = [kcum k];
        end
    end
    zothers = sum(z(:,1:nshocks,:),2);
    shock_groups.(['group' int2str(ngroups+1)]).label =  'Others';
    shock_groups.(['group' int2str(ngroups+1)]).shocks =  cellstr(M_.exo_names(find(~ismember([1:M_.exo_nbr],kcum)),:))';
    %M_.shock_groups.(options_.use_shock_groups)=shock_groups;
    if any(any(zothers)),
        shock_names = [shock_names; {'Others + Initial Values'}];
    end        
    zz(:,ngroups+1,:) = sum(z(:,1:nshocks+1,:),2);
    zz(:,ngroups+2,:) = z(:,nshocks+2,:);
    z = zz;
else
    shock_names = M_.exo_names;
end

%%

% Retrieves the variable names in TeX 
for j=1:size(varlist,1)
	iendo(j) =  strmatch(strtrim(varlist(j,:)),M_.endo_names,'exact');
    texvarlist{j}=deblank(M_.endo_names_tex(iendo(j),:));
end

% Retrieves the variable indices in M_ by their names
for j=1:size(varlist,1)
    indx(j) = find(strcmp(strtrim(varlist(j,:)),cellstr(M_.endo_names)));
end

% Retrieves the shock names
for i=1:ngroups
    for ii=1:length(M_.shock_groups.(options_.use_shock_groups).(sprintf('group%d', i)).shocks)
        shock_name_tmp = shock_groups.(sprintf('group%d', i)).shocks;
        indbuf = strmatch(shock_name_tmp{ii},M_.exo_names,'exact');
        if ~isempty(indbuf),
            index{i}(ii) = indbuf;
        elseif ~isempty(shock_names{i}{ii}),
            error(['Shock name ',shock_names{i}{ii}, ' not found.' ]);
        end
    end
end

% Lists the files of smoothed variables , innovations and parameters +
% stock_inno
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


% Set-up the time frame
% T: file listing the observation dates (or indicies)
load data T;
gend=size(stock,2);
if all(T==round(T))
    % annual model
    TT=0:gend;
else
    TT=[0:0.25:ceil(gend/4)];
end
if exist('T')
    TT=T(options_.first_obs)+TT;
end

%if nargin<5 || isempty(T0),
if(isempty(T0))
    t1decomp =25; % set it manually to change the initial point for plotting shock decomp
else
    switch type
        case {'qoq' , ''}
            t1decomp = max(1,find(TT==T0));
        case 'yoy'        
            t1decomp = max(5,find(TT==T0));
    end
end


switch type
    case {'qoq' , ''}
        kpoint = t1decomp:gend;
        gend2 = gend;
    case 'yoy'
        kpoint4 = t1decomp:gend;
        kpoint3 = kpoint4-1;
        kpoint2 = kpoint4-2;
        kpoint1 = kpoint4-3;
        gend2 = length(kpoint4);
end

% 
B=0;
SS=[];
if options_.debug
    save debug.mat
end
for j=1:length(param_file_list),
    load([DirectoryName,filesep,param_file_list(j).name],'stock_ys')
    stock_params=load([DirectoryName,filesep,param_file_list(j).name],'stock');
    stock_decomp=zeros(size(varlist,1),size(shock_names,1),gend2,size(stock_params.stock,1));
    B=B+size(stock_params.stock,1);
    SS=[SS;stock_ys(1:size(stock_params.stock,1),indx)];
    for i=1:size(stock_params.stock,1),
        jsmooth = jsmooth+1;
        if jsmooth>size(stock,3),
            jsmooth=1;
            jsmooth_file = jsmooth_file+1;
            load([DirectoryName,filesep,M_.fname,'_smooth',int2str(jsmooth_file),'.mat'],'stock')
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

        a1=ahat(:,1);
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
        inn=inn(oo_.dr.inv_order_var,:);
        deco=deco(oo_.dr.inv_order_var,:,:);
      
        switch type
            case {'qoq' , ''}  
                inno = inn(indx,kpoint);
                sdec0 = deco(indx,:,kpoint);                
            case 'yoy'     
                inno1 = inn(indx,kpoint1);
                inno2 = inn(indx,kpoint2);
                inno3 = inn(indx,kpoint3);
                inno4 = inn(indx,kpoint4);

                sdec1 = (deco(indx,:,kpoint1));
                sdec2 = (deco(indx,:,kpoint2));
                sdec3 = (deco(indx,:,kpoint3));
                sdec4 = (deco(indx,:,kpoint4));
                sdec0 = sdec1 + sdec2 + sdec3 + sdec4;
        end
        for jx=1:size(varlist,1)
            for ix=1:size(shock_names,1)-1
                sdec(jx,ix,:)=sum(sdec0(jx,index{ix},:),2);
                sdec0(jx,index{ix},:)=0;
            end
            switch type
                case {'qoq' , ''}
                    sdec_f(jx,:,:)=[squeeze(sdec(jx,:,:)); as(indx(jx),kpoint)-inno(jx,:)+sum(squeeze(sdec0(jx,:,:)),1)];                    
                case 'yoy'
                    dump2=as(indx(jx),kpoint1)+as(indx(jx),kpoint2)+as(indx(jx),kpoint3)+as(indx(jx),kpoint4)-(inno1(jx,:)+inno2(jx,:)+inno3(jx,:)+inno4(jx,:));
                    sdec_f(jx,:,:)=[squeeze(sdec(jx,:,:)); dump2+sum(squeeze(sdec0(jx,:,:)),1)];
            end
        end
        stock_decomp(:,:,:,i) = sdec_f;
    end
    stock0=stock;
    stock=stock_decomp;
    save([DirectoryName,filesep,M_.fname,'_decomp_',type,'_vs_',num2str(size(shock_names,1)),'Shocks',file_name,'_',int2str(j)],'stock')
    stock=stock0;
    clear stock0;
end

stock1 = zeros(size(varlist,1),size(shock_names,1),gend2,B);

k = 0;
for file = 1:length(param_file_list)
    load([DirectoryName,filesep,M_.fname,'_decomp_',type,'_vs_',num2str(size(shock_names,1)),'Shocks',file_name,'_',int2str(file)],'stock')
    k = k(end)+(1:size(stock,4));
    stock1(:,:,:,k) = stock;
end
clear stock    
SSMean = mean(SS);
for j = 1:size(varlist,1),
    stock2=sum(stock1,2);
    for is = 1:gend2
        [YMean(j,is),YMedian(j,is),YVar(j,is),YHPD(:,j,is),YDistrib(:,j,is)] = ...
            posterior_moments(squeeze(stock2(j,1,is,:)),0,options_.mh_conf_sig);
    end

    for i=1:size(shock_names,1),
        for is = 1:gend2
            [Mean(j,i,is),Median(j,i,is),Var(j,i,is),HPD(:,j,i,is),Distrib(:,j,i,is)] = ...
                posterior_moments(squeeze(stock1(j,i,is,:)),0,options_.mh_conf_sig);
        end
    end
end

ngroups = size(shock_names,1)+1;
nrow=round(sqrt(ngroups));
if nrow^2 < ngroups,
    ncol=nrow+1;
else
    ncol=nrow;
end


%% Plotting
for j=1:size(varlist,1),
    nbcmpts=size(Mean,2);
       
    if flagmap==0
        func = @(x) colorspace('RGB->Lab',x);
        MAP = distinguishable_colors(nbcmpts,'w',func);
        MAP(end,:) = [0.7 0.7 0.7];
    else
        colours=leg(:,2);
        
        for i=1:length(colours); 
            MAP(i,:)=rgb(colours{i,:}); 
        end
    end
    
    
    %leg0=leg(:,1); % leg0 => shock_names
    
    sdec2=squeeze(Mean(j,:,:));
    as=sum(sdec2,1);
    ipos2=sdec2>0;
    ineg2=sdec2<0;


    hinit = dyn_figure(options_.nodisplay,'name',strtrim(varlist(j,:)), 'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','landscape','renderermode','auto');
    colormap(MAP);
    set(gcf,'position' ,[50 50 1100 850])
    hp=bar((sdec2.*ipos2)','stacked','EdgeColor',[0 0 0]);
    shading faceted;
    %   set(get(hp(end),'children'),'Facecolor',[1 1 1]);
    hold on, hn=bar((sdec2.*ineg2)','stacked');
    shading faceted;
    %   set(get(hn(end),'children'),'Facecolor',[1 1 1]);
    
    set(gca,'position',[0.05 0.4 0.9 0.55],'units','normalized')
    hleg=legend(shock_names,'interpreter','none','location','Best');
    shading faceted;
    set(hleg,'position',[0.35 0.05 0.3 0.3],'units','normalized')
    
    hleg2 = get(hleg,'children');
    child = get(hleg2(1:2:end),'children');
    for ichild = 1:length(child),
        set(child{ichild},'Edgecolor',[0 0 0]);
    end
    hold on, 
    plot(squeeze(YHPD(:,j,:))',':k')
    switch type
        case {'qoq' , ''}
            h1=plot(as(kpoint),'k-d');
            tit2 = ['quarter to quarter ',strtrim(varlist(j,:))];
            set(gca,'Xtick',[1:4:length(kpoint)])
            set(gca,'Xticklabel',TT(kpoint(1):4:end))
            set(gca,'xlim',[0 length(kpoint)+4])            
        case 'yoy'
            h1=plot(as,'k-d');
            tit2 = ['year on year ',strtrim(varlist(j,:))];
            set(gca,'Xtick',[1:4:length(kpoint1)+4])
            set(gca,'Xticklabel',TT(kpoint4(1):4:end))
            set(gca,'xlim',[0 length(kpoint1)+4])
    end
    set(h1,'MarkerFaceColor', 'k')
    title(tit2,'interpreter','none')
    a=axis;
    
    switch type
        case {'qoq' , ''}
            SSMean_f = SSMean;
        case 'yoy'
            SSMean_f = SSMean*4;
    end
    if a(3)+SSMean_f(j)<0 && a(4)+SSMean_f(j)>0,
        plot(a(1:2),SSMean_f(j)*[-1 -1],'k--','linewidth',1)
        ytick=get(gca,'ytick');
        ytick1=ytick-SSMean_f(j);
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
    set(gca,'yticklabel',num2str(get(gca,'ytick')'+SSMean_f(j),'%4.2g'));
    dyn_saveas(gcf,[OutDirectoryName,filesep,M_.fname,'_posterior_shock_dec_', type,'_',strtrim(varlist(j,:)),'_vs_',num2str(size(shock_names,1)),'Shocks',file_name,'_Init'],options_.nodisplay,options_.graph_format)

    hinit = dyn_figure(options_.nodisplay,'name',strtrim(varlist(j,:)), 'PaperPositionMode', 'auto','PaperType','A4','PaperOrientation','landscape','renderermode','auto');
    
    switch type
        case {'qoq' , ''}
            for i=1:size(shock_names,1)
                distr=squeeze(Distrib(:,j,i,:));
                subplot(nrow,ncol,i), plot(distr','color',MAP(i,:)), title(shock_names{i}),
                set(gca,'Xtick',[1:16:length(kpoint)])
                set(gca,'Xticklabel',TT(kpoint(1):16:end))
                set(gca,'xlim',[0 length(kpoint)+4])
                hold all, plot(sdec2(i,:),'k')
            end
        case 'yoy'
            a0=a*0;
            a0(3)=inf;
            a0(4)=-inf;

            for i=1:size(shock_names,1)
                distr=squeeze(Distrib(:,j,i,:));
                subplot(nrow,ncol,i), 
                plot(distr','color',MAP(i,:)), title(shock_names{i}),
                hold all, plot(sdec2(i,:),'k'), 
                plot([0 length(kpoint1)+4],[0 0],'k--')
                plot(as, 'k', 'linewidth',2)
                axis tight;
                a=axis;        
                set(gca,'Xtick',[1:16:length(kpoint1)+4])
                set(gca,'Xticklabel',TT(kpoint4(1):16:end))
                set(gca,'xlim',[0 length(kpoint1)+4])
                a0(3)=min(a(3),a0(3));
                a0(4)=max(a(4),a0(4));
                set(gca,'ylim',a0(3:4))
            end
    end    
    dyn_saveas(gcf,[OutDirectoryName,filesep,M_.fname,'_posterior_shock_dec_', type,'_',strtrim(varlist(j,:)),'_vs_',num2str(size(shock_names,1)),'Shocks',file_name,'_Init_Detail'],options_.nodisplay,options_.graph_format)
    
    
    % Including the figure in the TeX document
    if options_.TeX
        
        tempfilename=[OutDirectoryName,filesep,M_.fname,'_posterior_shock_dec_',type,'_',strtrim(varlist(j,:)),'_vs_',num2str(size(shock_names,1)),'Shocks',file_name,'_Init'];
        
        [a,b,c]=fileparts(tempfilename);
        
        
        % TeX eps loader file
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.8][0]{%s}\n',tit2,['quarter on quarter $' deblank(texvarlist{j}) '$']);
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,['\\includegraphics[width=0.65\\textwidth,angle=-90] {' b '} \n']);
        fprintf(fidTeX,'\\caption{Shock Decomposition quarter on quarter $%s$}',texvarlist{j});
        fprintf(fidTeX,'\\label{Fig:Shock Decomposition:%s}\n',int2str(j));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
        
    end
    
end

if options_.TeX
    fclose(fidTeX);
end