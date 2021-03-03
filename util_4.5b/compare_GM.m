function compare_GM(folders,var_list,comp_type,varargin)
% HELP: This function compares different model specifications,
% which refer to different folders in the first input argument, along four
% dimensions which can be entered in the 3rd input argument, comp_type.
%
% You can enter either 'gemc_fit',or 'frcst', 'irfs'or 'post_dens'
%
% The second input argument contains the variables/parameters you are interested in comparing.
% If you choose 'gemc_fit' the sintax is going to be:
% compare_GM({'path1',..'pathn'},{'var1',..'var9'},'gemc_fit')
%
% If you choose 'frcst' the sintax is going to be:
% compare_GM({'path1',..'pathn'},{'var1',..'var9'},'frcst',nfrcst), where
% nfrcst is the number of quarters out of sample
%
% If you choose 'irfs' the sintax is going to be:
% compare_GM({'path1',..'pathn'},{'var1',..'varn'},'irfs','shockname'), where
% shockname  is the exogenous shock you are interested in.
%
% If you want to compare 'irfs' across countries:
% compare_GM({'path1',..'pathn'},{'var1_',..'varn_'},'irfs','shockname_',{co1,..,con}), where
% shockname  is the exogenous shock you are interested in and co1.. are the
% country codes
%
%
%If you choose 'post_dens' the sintax is going to be:
% compare_GM({'path1',..'pathn'},{'param1',..'paramn'},'post_dens'), 
% if you want to plot posterior densisties of ALL parameters you should
% invoke 
%compare_GM({'path1',..'pathn'},'all','post_dens'),
%if you want to compare post densities of parameters across countries the syntax is going to be
%compare_GM({'path1',..'pathn'},{'param1_',..'paramn_'},'post_dens',{'co1[_cox]',..,'con[_cox]'),
    


current_folder=pwd;

if ~isempty(strmatch('mode', comp_type, 'exact')) && length(folders)==1
    cd(folders{1})
    load gemc_mh_mode xparam1 parameter_names
    xparam1_1 = xparam1;
    parameter_names_1 = parameter_names; 
    load gemc_mh_mode0 xparam1 parameter_names
    xparam1_2 = xparam1;
    parameter_names_2 = parameter_names; 
    clear xparam1 parameter_names
    
elseif ~isempty(strmatch('mode', comp_type, 'exact')) && length(folders)==2
   for i=1:length(folders)
       cd(folders{i})
       tmp=load('gemc_mh_mode','xparam1','parameter_names');
       eval(['xparam1_',int2str(i),'= tmp.xparam1;']);
       eval(['parameter_names_',int2str(i),'= tmp.parameter_names;']);
       cd(current_folder)
   end
       [parameter_names_1,idx1]=sort(parameter_names_1);
       [parameter_names_2,idx2]=sort(parameter_names_2);
       xparam1_1 = xparam1_1(idx1);
       xparam1_2 = xparam1_2(idx2);
       length1 = length(parameter_names_1);
       length2 = length(parameter_names_2);
       if length1 ~= length2
           if length1>length2
               for j=length2+1:length1
              parameter_names_2{j}='missing';
              xparam1_2(j)=NaN;
               end
           else
               for j=length1+1:length2
              parameter_names_1{j}='missing';
              xparam1_1(j)=NaN; 
               end
           end
       end
       [common1, icommon1] = ismember(parameter_names_1,parameter_names_2);
       [common2, icommon2] = ismember(parameter_names_2,parameter_names_1);
       parameter_names_1=parameter_names_1(common1);
       xparam1_1=xparam1_1(common1);
       parameter_names_2=parameter_names_2(common2);
       xparam1_2=xparam1_2(common2);
end
if ~isempty(strmatch('mode', comp_type, 'exact')) && length(folders)>2
     error(sprintf('Maximum number of folders is 2 when comparing modes'))    
end
% if strmatch('mode', comp_type, 'exact')==0     
if isempty(strmatch('mode', comp_type, 'exact'))    
for i=1:length(folders)
    cd(folders{i})
    if strmatch('frcst', comp_type, 'exact')
        tmp0=load('gemc_results_cfrcst','oo_');
        tmp=load('gemc_results','M_','options_');
        tmp.oo_=tmp0.oo_;
        clear tmp0;
    else
        if strmatch('data', comp_type, 'exact')
          tmp=load('gemc_results','options_'); 
          tmp.dataobs= load('dataobs');
        else
        tmp=load('gemc_results','oo_','M_','options_');
        end
    end
    indx = strfind(folders{i},'_');
    if length(indx)>6
    lgnd{i}= folders{i}(indx(end-6)+1:end);
    else
    lgnd{i}= folders{i};
    end
    if strmatch('data', comp_type, 'exact')
        eval(['options',int2str(i),'_=tmp.options_;']);
        eval(['dataobs',int2str(i),'=tmp.dataobs;']);
    else
    eval(['oo',int2str(i),'_= tmp.oo_;']);
    eval(['M',int2str(i),'_=tmp.M_;']);
    eval(['options',int2str(i),'_=tmp.options_;']);
    end
    %eval(['bayestopt',int2str(i),'_=bayestopt_;']);
    options_ = tmp.options_;
    clear tmp
    cd(current_folder)
end


if exist(options1_.datafile)
    instr = options1_.datafile;
else
    instr = ['load ' options1_.datafile];
end
try
    eval([instr ' T']);
    if ~exist('T','var'),
        temp = eval(deblank(options1_.varobs{1}));
        T=[1:length(temp)]';
        clear temp;
    end
catch
    %T=[1:options1_.nobs-options1_.forecast];
    T=[1:options1_.first_obs+options1_.nobs+options1_.forecast-1];
end
end
if  strmatch('gemc_smooth_var', comp_type, 'exact')
    
    % options1_.nobs=options1_.nobs-options1_.forecast;
    
    fobs = options1_.first_obs;
    nfig=ceil(length(var_list)/9);
    for ifig=1:nfig
        if ifig<nfig
            nplots=9;
        else
            nplots= rem(length(var_list),9);
            if nplots==0
                nplots=9;
            end
        end
        figure(ifig)
    
    for j=1:nplots
        subplot(3,3,j)
        for i=1:length(folders)
            eval(['plot(T(fobs:fobs+options1_.nobs-1),oo',int2str(i),'_.SmoothedVariables.',var_list{j+(ifig-1)*9},')'])
            hold on
        end
        
        hh=get(gca,'children');
        ym = inf;
        yM = -inf;
        for ih=1:length(hh),
            try
                ym = min(ym,min(get(hh(ih),'ydata')));
                yM = max(yM,max(get(hh(ih),'ydata')));
            catch
            end
        end
        
        %        hold on, plot([T(fobs)-1 T(fobs+options1_.nobs-1)+1],[get_mean(var_list(1,:)) get_mean(var_list(1,:))],'r--')
        set(gca,'xlim',[T(fobs)-1 T(fobs+options1_.nobs-1)+1])
        dylim=max((yM-ym)*0.05,1.e-10);
        set(gca,'ylim',[ym-dylim yM+dylim])
        title(var_list{j+(ifig-1)*9},'interpreter','none')
    end
    end
    legend(lgnd,'interpreter','none')
end

if  strmatch('gemc_smooth_shocks', comp_type, 'exact')
    
    % options1_.nobs=options1_.nobs-options1_.forecast;
    
    fobs = options1_.first_obs;
    nfig=ceil(length(var_list)/9);
    for ifig=1:nfig
        if ifig<nfig
            nplots=9;
        else
            nplots= rem(length(var_list),9);
            if nplots==0
                nplots=9;
            end
        end
        figure(ifig)
    
    for j=1:nplots
        subplot(3,3,j)
        for i=1:length(folders)
            eval(['plot(T(end-11:end-4),oo',int2str(i),'_.SmoothedShocks.',var_list{j+(ifig-1)*9},'(end-7:end)',')'])
            hold on
        end
        
        hh=get(gca,'children');
        ym = inf;
        yM = -inf;
        for ih=1:length(hh),
            try
                ym = min(ym,min(get(hh(ih),'ydata')));
                yM = max(yM,max(get(hh(ih),'ydata')));
            catch
            end
        end
        
        %        hold on, plot([T(fobs)-1 T(fobs+options1_.nobs-1)+1],[get_mean(var_list(1,:)) get_mean(var_list(1,:))],'r--')
        set(gca,'xlim',[T(end-11)-1 T(end-4)+1])
        dylim=max((yM-ym)*0.05,1.e-10);
        set(gca,'ylim',[ym-dylim yM+dylim])
        title(var_list{j+(ifig-1)*9},'interpreter','none')
    end
    end
    legend(lgnd,'interpreter','none')
end

if   strmatch('gemc_fit', comp_type, 'exact')
    
    % options1_.nobs=options1_.nobs-options1_.forecast;
    
    istart = max(2,options1_.presample+1);
    fobs = options1_.first_obs;
    nplots= length(var_list);
    for j=1:nplots
        subplot(3,3,j)
        for i=1:length(folders)
            eval(['plot(T(fobs+istart-1:fobs+options1_.nobs-1),oo',int2str(i),'_.FilteredVariables.',var_list{j},'(istart-1:options1_.nobs-1) )'])
            hold on
        end
        
        eval(['plot(T(fobs+istart-1:fobs+options1_.nobs-1),oo1_.SmoothedVariables.',var_list{j},'(istart:end),''color'',[0.7 0.7 0.7])'])
        %eval(['plot(T(fobs+istart-1:fobs+options1_.nobs-1),oo1_.UpdatedVariables.',var_list{j},'(istart:end),''color'',[0.7 0.7 0.7])'])
        hh=get(gca,'children');
        ym = inf;
        yM = -inf;
        for ih=1:length(hh),
            try
                ym = min(ym,min(get(hh(ih),'ydata')));
                yM = max(yM,max(get(hh(ih),'ydata')));
            catch
            end
        end
        
        %        hold on, plot([T(fobs)-1 T(fobs+options1_.nobs-1)+1],[get_mean(var_list(1,:)) get_mean(var_list(1,:))],'r--')
        set(gca,'xlim',[T(fobs)-1 T(fobs+options1_.nobs-1)+1])
        dylim=max((yM-ym)*0.05,1.e-10);
        set(gca,'ylim',[ym-dylim yM+dylim])
        title(var_list{j},'interpreter','none')
    end
    legend(lgnd,'interpreter','none')
    for i=1:length(folders)
        %       oo_=eval(['oo',int2str(j),'_;']);
        %       M_=eval(['M',int2str(j),'_;']);
        %       options_=eval(['options',int2str(j),'_;']);
        %       bayestopt_=eval(['bayestopt',int2str(j),'_;']);
        %
        %      [rmse, lnam, r2(j,:) ,vv] = plot_fit(var_list{:});
        lgy_ = eval(['M',int2str(i),'_.endo_names;']);
        mfys=[];
        for j=1:length(var_list),
            dum = strmatch(var_list{j},lgy_,'exact');
            mfys = [mfys dum];
            if j==1,
                lgobs_ = var_list{j};
            else
                lgobs_ = str2mat(lgobs_,var_list{j});
            end
        end
        ys_ = eval(['oo',int2str(i),'_.steady_state;']);
        if eval(['options',int2str(i),'_.loglinear == 1;'])
            constant = log(ys_(mfys));
        else
            constant = ys_(mfys);
        end
        trend = eval(['constant*ones(1,options1_.nobs);']);
        
        %         trend = eval(['constant*ones(1,options',int2str(i),'_.nobs);']);
        for j=1:length(var_list)
            eval(['vv(:,j) = (oo',int2str(i),'_.UpdatedVariables.',var_list{j},'(istart:end)-oo',int2str(i),'_.FilteredVariables.',var_list{j},'(istart-1:options1_.nobs-1));'])
            eval(['r2(j,i) = 1-sum(vv(:,j).^2)/sum((oo',int2str(i),'_.UpdatedVariables.',var_list{j},'(istart:end)-trend(j,istart:end)'').^2);'])
        end
    end
    
    T = table(r2,'RowNames',var_list');
    save R2 T
end

if strmatch('frcst', comp_type, 'exact')
    if nargin==4
        options1_.forecast=varargin{1};
    else
        options1_.forecast=8;
    end
    %options1_.nobs = options1_.nobs + options1_.forecast;
    if eval(['~isfield(oo',int2str(i),'_, ''q2avec'')'])
        co={};
        for j=1:length(var_list)
            indx = strfind(var_list{j},'_');
            co_tmp= var_list{j}(indx(end)+1:end);
            if j==1
                co{1}=co_tmp;
                z=1;
            else
                if strmatch(co_tmp, co,'exact')
                else
                    z=z+1;
                    co{z}=co_tmp;
                end
                
            end
            
        end
        avname=0;
            q2avec=0;
        for j=1:z
            M_=M1_;
            co_tmp=co{j};
            
            if isint(avname)
                avname = strcat('LYOBS_',co_tmp);
            else
                avname = char(avname,strcat('LYOBS_',co_tmp));
            end
            q2a=struct();
            q2a.type=1;
            q2a.islog=1;
            i = strmatch('GYTREND0',M_.param_names,'exact');
            
            if isempty(i)
                error(sprintf('Can''t find parameter GYTREND0'))
            end
            q2a.GYTREND0 = M_.params(i);
            q2a.cumfix = 1;
            q2a.plot = 1;
            q2a.aux = 0;
            q2a.name=strcat('LYOBSA_',co_tmp);
            q2a.gname=strcat('GYOBSA_',co_tmp);
            if isint(q2avec)
                q2avec=q2a;
            else
                q2avec=[q2avec q2a];
            end
            avname = char(avname,strcat('LPYOBS_',co_tmp));
            q2a=struct();
            q2a.type=2;
            q2a.islog=1;
            q2a.aux = 0;
            i = strmatch('GP0',M_.param_names,'exact');
            
            if isempty(i)
                error(sprintf('Can''t find parameter %s', pname))
            end
            q2a.GYTREND0 = M_.params(i);
            q2a.cumfix = 1;
            q2a.plot = 1;
            q2a.name=strcat('LPYOBSA_',co_tmp);
            q2a.gname=strcat('PHIYOBSA_',co_tmp);
            q2avec=[q2avec q2a];
            avname = char(avname,strcat('INOM_',co_tmp));
            q2a=struct();
            q2a.type=1;
            q2a.islog=0;
            q2a.aux = 0;
            q2a.GYTREND0 = 0;
            q2a.cumfix = 1;
            q2a.plot = 2;
            q2a.name=strcat('INOMA_',co_tmp);
            q2a.gname='';
            q2avec=[q2avec q2a];
            avname = char(avname,strcat('LE_',co_tmp));
            q2a=struct();
            q2a.type=2;
            q2a.islog=1;
            q2a.aux = 0;
            q2a.GYTREND0 = 0;
            q2a.cumfix = 1;
            q2a.plot = 1;
            q2a.name=strcat('LEA_',co_tmp);
            q2a.gname=strcat('GEA_',co_tmp);
            q2avec=[q2avec q2a];
            if ~strcmp(co_tmp,'RoW') && ~strcmp(co_tmp,'REA')
                avname = char(avname,strcat('LI_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=1;
                i = strmatch('GYTREND0',M_.param_names,'exact');
                ii= strmatch(strcat('GAPI0_',co_tmp),M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    if min([~strcmp('DE',co) ~strcmp('IT',co) ~strcmp('ES',co) ~strcmp('FR',co)])==1 %|| ~strcmp('IT',co) || ~strcmp('ES',co) || ~strcmp('FR',co)
                    error(sprintf('Can''t find parameter  GYTREND0 or GAPI0'))
                    end
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('LIA_',co_tmp);
                q2a.gname=strcat('GIA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LPI_',co_tmp));
                q2a=struct();
                q2a.type=2;
                q2a.islog=1;
                q2a.aux = 0;
                i = strmatch('GP0',M_.param_names,'exact');
                ii= strmatch(strcat('GAPI0_',co_tmp),M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    if min([~strcmp('DE',co) ~strcmp('IT',co) ~strcmp('ES',co) ~strcmp('FR',co)])==1 %|| ~strcmp('IT',co) || ~strcmp('ES',co) || ~strcmp('FR',co)
                    
                    error(sprintf('Can''t find parameter  GP0 or GAPI0'))
                    end
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.name=strcat('LPIA_',co_tmp);
                q2a.gname=strcat('PHIIA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LIG_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=1;
                i = strmatch('GYTREND0',M_.param_names,'exact');
                ii= strmatch(strcat('GAPI0_',co_tmp),M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    if min([~strcmp('DE',co) ~strcmp('IT',co) ~strcmp('ES',co) ~strcmp('FR',co)])==1 %|| ~strcmp('IT',co) || ~strcmp('ES',co) || ~strcmp('FR',co)
                    
                    error(sprintf('Can''t find parameter  GYTREND0 or GAPI0'))
                    end
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('LIGA_',co_tmp);
                q2a.gname=strcat('GIGA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LPIG_',co_tmp));
                q2a=struct();
                q2a.type=2;
                q2a.islog=1;
                q2a.aux = 0;
                i = strmatch('GP0',M_.param_names,'exact');
                ii= strmatch(strcat('GAPI0_',co_tmp),M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    if min([~strcmp('DE',co) ~strcmp('IT',co) ~strcmp('ES',co) ~strcmp('FR',co)])==1 %|| ~strcmp('IT',co) || ~strcmp('ES',co) || ~strcmp('FR',co)
                    
                    error(sprintf('Can''t find parameter  GYTREND0 or GAPI0'))
                    end
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.name=strcat('LPIGA_',co_tmp);
                q2a.gname=strcat('PHIIGA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LC_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=1;
                i = strmatch('GYTREND0',M_.param_names,'exact');
                ii= strmatch(strcat('GAPC0_',co_tmp),M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    if min([~strcmp('DE',co) ~strcmp('IT',co) ~strcmp('ES',co) ~strcmp('FR',co)])==1 %|| ~strcmp('IT',co) || ~strcmp('ES',co) || ~strcmp('FR',co)
                    
                    error(sprintf('Can''t find parameter  GYTREND0 or GAPC0'))
                    end
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('LCA_',co_tmp);
                q2a.gname=strcat('GCA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LPCVAT_',co_tmp));
                q2a=struct();
                q2a.type=2;
                q2a.islog=1;
                q2a.aux = 0;
                i = strmatch('GP0',M_.param_names,'exact');
                ii= strmatch(strcat('GAPC0_',co_tmp),M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    if min([~strcmp('DE',co) ~strcmp('IT',co) ~strcmp('ES',co) ~strcmp('FR',co)])==1 %|| ~strcmp('IT',co) || ~strcmp('ES',co) || ~strcmp('FR',co)
                    
                    error(sprintf('Can''t find parameter  GP0 or GAPC0'))
                    end
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.name=strcat('LPCVATA_',co_tmp);
                q2a.gname=strcat('PHICVATA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LN_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=1;
                i = strmatch('GPOP0',M_.param_names,'exact');
                
                if isempty(i)
                    error(sprintf('Can''t find parameter GPOP0'))
                end
                q2a.GYTREND0 = M_.params(i);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('LNA_',co_tmp);
                q2a.gname=strcat('GNA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LW_',co_tmp));
                q2a=struct();
                q2a.type=2;
                q2a.islog=1;
                q2a.aux = 0;
                i = strmatch('GP0',M_.param_names,'exact');
                ii= strmatch('GYTREND0',M_.param_names,'exact');
                iii= strmatch('GPOP0',M_.param_names,'exact');
                if isempty(i) || isempty(ii) || isempty(iii)
                    error(sprintf('Can''t find parameter  GPOP0 or GYTREND0 or GPOP0'))
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii)-M_.params(iii);
                
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.name=strcat('LWA_',co_tmp);
                q2a.gname=strcat('PHIWA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LG_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=1;
                i = strmatch('GYTREND0',M_.param_names,'exact');
                ii= strmatch(strcat('GAPG0_',co_tmp),M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    if min([~strcmp('DE',co) ~strcmp('IT',co) ~strcmp('ES',co) ~strcmp('FR',co)])==1 %|| ~strcmp('IT',co) || ~strcmp('ES',co) || ~strcmp('FR',co)
                    
                    error(sprintf('Can''t find parameter  GYTREND0 or GAPG0'))
                    end
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('LGA_',co_tmp);
                q2a.gname=strcat('GCGA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LPG_',co_tmp));
                q2a=struct();
                q2a.type=2;
                q2a.islog=1;
                q2a.aux = 0;
                i = strmatch('GP0',M_.param_names,'exact');
                ii= strmatch(strcat('GAPG0_',co_tmp),M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    if min([~strcmp('DE',co) ~strcmp('IT',co) ~strcmp('ES',co) ~strcmp('FR',co)])==1 %|| ~strcmp('IT',co) || ~strcmp('ES',co) || ~strcmp('FR',co)
                    
                    error(sprintf('Can''t find parameter  GP0 or GAPG0'))
                    end
                end
                q2a.GYTREND0 = M_.params(i)-M_.params(ii);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.name=strcat('LPG_',co_tmp);
                q2a.gname=strcat('PHIGA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('MTOT_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=0;
                i = strmatch('GYTREND0',M_.param_names,'exact');
                
                if isempty(i)
                    error(sprintf('Can''t find parameter  GYTREND0 '))
                end
                q2a.GYTREND0 = M_.params(i);
                
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('MTOTA_',co_tmp);
                q2a.gname=strcat('GMTOTA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LPMTOT_',co_tmp));
                q2a=struct();
                q2a.type=2;
                q2a.islog=1;
                q2a.aux = 0;
                i = strmatch('GP0',M_.param_names,'exact');
                
                if isempty(i)
                    error(sprintf('Can''t find parameter  GP0 '))
                end
                q2a.GYTREND0 = M_.params(i);
                
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.name=strcat('LPMTOTA_',co_tmp);
                q2a.gname=strcat('PHIMTOTA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LX_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=1;
                i = strmatch('GYTREND0',M_.param_names,'exact');
                
                if isempty(i)
                    error(sprintf('Can''t find parameter  GYTREND0 '))
                end
                q2a.GYTREND0 = M_.params(i);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('LXA_',co_tmp);
                q2a.gname=strcat('GXA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('LPX_',co_tmp));
                q2a=struct();
                q2a.type=2;
                q2a.islog=1;
                q2a.aux = 0;
                i = strmatch('GP0',M_.param_names,'exact');
                
                if isempty(i)
                    error(sprintf('Can''t find parameter  GP0 '))
                end
                q2a.GYTREND0 = M_.params(i);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.name=strcat('LPXA_',co_tmp);
                q2a.gname=strcat('PHIXA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('GMTOTN_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=2;
                i = strmatch('GYTREND0',M_.param_names,'exact');
                ii= strmatch('GP0',M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    error(sprintf('Can''t find parameter  GYTREND0 or GP0'))
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('LMTOTNA_',co_tmp);
                q2a.gname=strcat('GMTOTNA_',co_tmp);
                q2avec=[q2avec q2a];
                avname = char(avname,strcat('GXN_',co_tmp));
                q2a=struct();
                q2a.type=1;
                q2a.islog=2;
                i = strmatch('GYTREND0',M_.param_names,'exact');
                ii= strmatch('GP0',M_.param_names,'exact');
                if isempty(i) || isempty(ii)
                    error(sprintf('Can''t find parameter  GYTREND0 or GP0'))
                end
                q2a.GYTREND0 = M_.params(i)+M_.params(ii);
                q2a.cumfix = 1;
                q2a.plot = 1;
                q2a.aux = 0;
                q2a.name=strcat('LXNA_',co_tmp);
                q2a.gname=strcat('GXNA_',co_tmp);
                q2avec=[q2avec q2a];
            end
            avname = char(avname,strcat('TB_',co_tmp));
            q2a=struct();
            q2a.type=1;
            q2a.islog=0;
            i = strmatch('GYTREND0',M_.param_names,'exact');
            ii= strmatch('GP0',M_.param_names,'exact');
            if isempty(i) || isempty(ii)
                error(sprintf('Can''t find parameter  GYTREND0 or GP0'))
            end
            q2a.GYTREND0 = M_.params(i)+M_.params(ii);
            q2a.cumfix = 1;
            q2a.plot = 2;
            q2a.aux=0;
            q2a.name=strcat('TBA_',co_tmp);
            q2a.gname='';
            q2avec=[q2avec q2a];
            avname = char(avname,strcat('NFA_',co_tmp));
            q2a=struct();
            q2a.type=1;
            q2a.islog=0;
            i = strmatch('GYTREND0',M_.param_names,'exact');
            ii= strmatch('GP0',M_.param_names,'exact');
            if isempty(i) || isempty(ii)
                error(sprintf('Can''t find parameter  GYTREND0 or GP0'))
            end
            q2a.GYTREND0 = M_.params(i)+M_.params(ii);
            q2a.cumfix = 1;
            q2a.plot = 2;
            q2a.aux=0;
            q2a.name=strcat('NFAA_',co_tmp);
            q2a.gname='';
            q2avec=[q2avec q2a];
        end
    end
    q2a=q2avec;
    
    var_name = cellstr(avname);
    T = T(options1_.first_obs:end);
    
    t0=min(find(T==floor(T)));
    forecast_ = ceil(options1_.forecast/4);
    T = T(t0:4:end);
    Tplot = T;
    TF  = T(end - forecast_ +1 : end);
    nplots= length(var_list);
    Tlim=[2005 2020];
    is     = find(T==Tlim(1));
    npreamble = (forecast_+1)*4-options1_.forecast;
    i = 1;
    for j=1:nplots
        oo_=eval(['oo',int2str(i),'_']);
        M_=eval(['M',int2str(i),'_']);
        options_=eval(['options',int2str(i),'_']);
        z = strmatch(var_list{j}, var_name);
        if isempty(z)
            error(sprintf('You can enter only these variables: %s', avname'))
            
        end
        [ya, yass, gya, gyass] = ...
            eval(['quarterly2annual(oo',int2str(i),'_.SmoothedVariables.',var_list{j},'(t0:end)- oo_.dr.ys(strmatch(''',var_list{j},''',M_.endo_names,''exact'')) , oo_.dr.ys(strmatch(''',var_list{j},''',M_.endo_names,''exact'')) ,q2a(',int2str(z),').GYTREND0,q2a(',int2str(z),').type,q2a(',int2str(z),').islog,q2a(',int2str(z),').aux);']);
        if q2a(z).plot ==1,
            ztmp=gya+gyass;
        elseif q2a(z).plot == 2
            ztmp=ya+yass;
        else
            error('choose level or growth rate')
        end
        EC_plot            = ztmp;
        actual_values      = EC_plot(1:length(T) - forecast_);
        steadyvar(j)=gyass;
        TimeLineQA = [is:length(EC_plot)];
        subplot(3,3,j)
        
        for i=1:length(folders)
           
            
            
            oo_=eval(['oo',int2str(i),'_']);
            M_=eval(['M',int2str(i),'_']);
            options_=eval(['options',int2str(i),'_']);
            ypreamble = eval(['oo',int2str(i),'_.SmoothedVariables.',var_list{j},'(end-(forecast_+1)*4+1:end);']);
            ypreamble = ypreamble(1:npreamble);
            % unconditional forecasts
            ytmp = eval(['[ypreamble; oo',int2str(i),'_.forecast.Mean.',var_list{j},']- oo_.dr.ys(strmatch(''',var_list{j},''',M_.endo_names,''exact'')) ;']);
            [ya, yass, gya, gyass] = ...
                quarterly2annual(ytmp, oo_.dr.ys(strmatch(var_list{j},M_.endo_names,'exact')) ,q2a(z).GYTREND0,q2a(z).type,q2a(z).islog,q2a(z).aux);
            
            %     [ya, yass, gya, gyass] = ...
            %         eval(['quarterly2annual(ytmp,get_mean(''',var_list{j},'''),q2a(',int2str(z),').GYTREND0,q2a(',int2str(z),').type,q2a(',int2str(z),').islog,q2a(',int2str(z),').aux);']);
            
            if q2a(z).plot ==1,
                ztmp=gya(2:end)+gyass;
            elseif q2a(z).plot == 2
                ztmp=ya(2:end)+yass;
            else
                error('choose level or growth rate')
            end
            forecastedvar(:,j,1) = ztmp;
            %%% this is for confidence intervals but it is not yet working
            %%% properly
            %             ytmp = eval(['[ypreamble; oo',int2str(i),'_.forecast.HPDsup.',var_list{j},']-get_mean(''',var_list{j},''');']);
            %             [ya, yass, gya, gyass] = ...
            %                 quarterly2annual(ytmp,get_mean(var_list{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
            %
            %             if q2a(j).plot ==1,
            %                 ztmp=gya(2:end)+gyass;
            %             elseif q2a(j).plot == 2
            %                 ztmp=ya(2:end)+yass;
            %             else
            %                 error('choose level or growth rate')
            %             end
            %             forecastedvar(:,j,2) = ztmp;
            %
            %             ytmp =eval(['[ypreamble; oo',int2str(i),'_.forecast.HPDinf.',var_list{j},']-get_mean(''',var_list{j},''');']);
            %             [ya, yass, gya, gyass] = ...
            %                 quarterly2annual(ytmp,get_mean(var_list{j}),q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
            %             if q2a(j).plot ==1,
            %                 ztmp=gya(2:end)+gyass;
            %             elseif q2a(j).plot == 2
            %                 ztmp=ya(2:end)+yass;
            %             else
            %                 error('choose level or growth rate')
            %             end
            %             forecastedvar(:,j,3) = ztmp;
            
            vplot           =  [actual_values ; forecastedvar(:,j,1)];
            % DSGE uncoditional forecast, upper bound
            % vplotStd1       = [actual_values ;  forecastedvar(:,j,2)];
            % DSGE uncoditional forecast, lower bound
            % vplotStd_1      = [actual_values ;  forecastedvar(:,j,3)];
            
            %             h = area([Tplot(TimeLineQA)],[vplotStd_1(is:end), vplotStd1(is:end) - vplotStd_1(is:end)]);
            %             set(h(2),'FaceColor',[.95 .95 .95])
            %             set(h(1),'FaceColor',[1 1 1])
            %             hold on;
            %%% this is for mean forecast confidence interval and it is not
            %%% working properly
            %             if eval(['isfield(oo',int2str(i),'_,''MeanForecast'')'])
            %                 ypreamble = eval(['oo',int2str(i),'_.SmoothedVariables.',var_list{j},'(end-(forecast_+1)*4+1:end);']);
            %                 ypreamble = ypreamble(1:npreamble);
            %
            %                 [ya, yass, gya, gyass] = ...
            %                     eval(['quarterly2annual([ypreamble; oo',int2str(j),'_.MeanForecast.HPDinf.',var_list{j},']-get_mean(''',var_list{j},'''),get_mean(''',var_list{j},'''),q2a(',int2str(z),').GYTREND0,q2a(',int2str(z),').type,q2a(',int2str(z),').islog,q2a(',int2str(z),').aux);']);
            %                 if q2a(j).plot ==1,
            %                     zinf=gya(2:end)+gyass;
            %                 elseif q2a(j).plot == 2
            %                     zinf=ya(2:end)+yass;
            %                 else
            %                     error('choose level or growth rate')
            %                 end
            %                 [ya, yass, gya, gyass] = ...
            %                     eval(['quarterly2annual([ypreamble; oo',int2str(j),'_.MeanForecast.HPDinf.',var_list{j},']-get_mean(''',var_list{j},'''),get_mean(''',var_list{j},'''),q2a(',int2str(z),').GYTREND0,q2a(',int2str(z),').type,q2a(',int2str(z),').islog,q2a(',int2str(z),').aux);']);
            %                 if q2a(j).plot ==1,
            %                     zsup=gya(2:end)+gyass;
            %                 elseif q2a(j).plot == 2
            %                     zsup=ya(2:end)+yass;
            %                 else
            %                     error('choose level or growth rate')
            %                 end
            %                 fill([Tplot(TimeLineQA(end-forecast_+1:end)); Tplot(TimeLineQA(end:-1:end-forecast_+1))],[zinf;zsup(end:-1:1)],[0.75 0.75 0.95]);
            %
            %             end
            save debug.mat
            hold on, plot(Tplot(TimeLineQA),vplot(is:end),'Linewidth',2)
        end
         hold on, plot(Tplot(TimeLineQA),EC_plot(is:end),'k-.','Linewidth',2)
            hold on, plot([Tplot(is) Tplot(end)],[steadyvar(j) steadyvar(j)],'k:')
        if q2a(z).plot ==1,
            title(q2a(z).gname,'interpreter','none')
        elseif q2a(z).plot == 2
            title(q2a(z).name,'interpreter','none')
        else
            error('choose level or growth rate')
        end
        %     T=T(options1_.first_obs:end);
        %     is     = find(T==2011);
        %
        %
        %     for j=1:nplots
        %
        %          smoothedvar=eval(['oo',int2str(i),'_.SmoothedVariables.',var_list{j},';']);
        %          actual_values      = smoothedvar(1:length(T) - options1_.forecast);
        %          TimeLineQA = [is+3:4:length(actual_values)+ options1_.forecast]-3;
        %          subplot(3,3,j)
        %
        %
        %         for i=1:length(folders)
        %             forecastedvar=eval(['oo',int2str(i),'_.forecast.Mean.',var_list{j},';']);
        %             vplot           =  [actual_values ; forecastedvar];
        %             if findstr('QA',var_list{j})~=0
        %               plot(T(TimeLineQA),vplot(is+3:4:end),'Linewidth',2)
        %             else
        %                plot(T(is:end),vplot(is:end),'Linewidth',2)
        %             end
        %             hold on
        %         end
        %         if findstr('QA',var_list{j})~=0
        %                 hold on, plot(T(TimeLineQA),smoothedvar(is+3:4:end),'k-.','Linewidth',2)
        %             else
        %                 hold on, plot(T(is:end),smoothedvar(is:end),'k-.','Linewidth',2)
        %             end
        %         title(var_list{j},'interpreter','none')
        %     end
        
    end
    legend(lgnd,'interpreter','none')
    
end

if strmatch('post_irfs', comp_type, 'exact')
    if nargin>=4
        shockid=varargin{1};
    end
    if nargin>=6
        periods=varargin{3};
    else
    periods=40;
    end
    if nargin>=7
        color=varargin{4};
    else
        color='color';
    end
    if strcmp(color,'BW')
        set(groot,'defaultAxesColorOrder',[0 0 0])
        set(groot,'defaultAxesLineStyleOrder',{'-','--',':','-.','-o','-x','-s'})
    end
    
    
    nfig=ceil(length(var_list)/9);
    for ifig=1:nfig
        if ifig<nfig
            nplots=9;
        else
            nplots= rem(length(var_list),9);
            if nplots==0
                nplots=9;
            end
        end
        figure(ifig)
        for j=1:nplots
        subplot(3,3,j)
            for i=1:length(folders)
                tmp=isletter(varargin{1});
                if tmp(end)
                    IRF=eval(['oo',int2str(i),'_.PosteriorIRF.dsge.Mean.',char(strcat(var_list{j+(ifig-1)*9},varargin{2}(i))),'_',shockid,';']);
                else
                    IRF=eval(['oo',int2str(i),'_.PosteriorIRF.dsge.Mean.',char(strcat(var_list{j+(ifig-1)*9},varargin{2}(i))),'_',char(strcat(shockid,varargin{2}(i))),';']);
                end
            hold on
            plot(1:periods,IRF(1:periods),'Linewidth',2)
            title(var_list{j+(ifig-1)*9},'interpreter','none')
            end
        end
        legend(lgnd,'interpreter','none')
    end
    if strcmp(color,'BW')
        set(groot,'defaultAxesColorOrder','remove')
        set(groot,'defaultAxesLineStyleOrder','remove')
    end
end

if strmatch('irfs', comp_type, 'exact')
    if nargin>=4
        shockid=varargin{1};
    end
    if nargin>=6
        periods=varargin{3};
    else
        periods=40;
    end
    if nargin>=7
        color=varargin{4};
    else
        color='color';
    end
    if nargin>=8
        jrc=[varargin{5} '.'];
    else
        jrc='';
    end
    if strcmp(color,'BW')
        set(groot,'defaultAxesColorOrder',[0 0 0])
        set(groot,'defaultAxesLineStyleOrder',{'-','--',':','-.','-o','-x','-s'})
    end
    nfig=ceil(length(var_list)/9);
    
    %check expansionary
    tmp=isletter(varargin{1});
    if tmp(end)
        tmp2=isletter(var_list{1});
        if tmp2(end)
            isexpansionary=eval(['oo',int2str(1),'_.',jrc,'irfs.',char(strcat('LYOBS_',var_list{1}(end-1:end))),'_',shockid,';']);
        else
            isexpansionary=eval(['oo',int2str(1),'_.',jrc,'irfs.',char(strcat('LYOBS_',varargin{2}(1))),'_',shockid,';']);
        end
    else
        tmp2=isletter(var_list{1});
        if tmp2(end)
            isexpansionary=eval(['oo',int2str(1),'_.',jrc,'irfs.',char(strcat('LYOBS_',var_list{1}(end-1:end))),'_',char(strcat(shockid,varargin{2}(1))),';']);
        else
            isexpansionary=eval(['oo',int2str(1),'_.',jrc,'irfs.',char(strcat('LYOBS_',varargin{2}(1))),'_',char(strcat(shockid,varargin{2}(1))),';']);
        end
        
    end
    isexpansionary=cumsum(isexpansionary(1:20));
    if isexpansionary(end)<0
       flip=-1;
    else
        flip=1;
    end
    for ifig=1:nfig
        if ifig<nfig
            nplots=9;
        else
            nplots= rem(length(var_list),9);
            if nplots==0
                nplots=9;
            end
        end
        figure(ifig)
        for j=1:nplots
        subplot(3,3,j)
        
            for i=1:length(folders)
                tmp=isletter(varargin{1});
                if tmp(end)
                    tmp2=isletter(var_list{j+(ifig-1)*9});
                    if tmp2(end)
                        IRF=eval(['oo',int2str(i),'_.',jrc,'irfs.',char(var_list{j+(ifig-1)*9}),'_',shockid,';']);
                    else
                        IRF=eval(['oo',int2str(i),'_.',jrc,'irfs.',char(strcat(var_list{j+(ifig-1)*9},varargin{2}(i))),'_',shockid,';']);
                    end
                else
                    tmp2=isletter(var_list{j+(ifig-1)*9});
                    if tmp2(end)
                        IRF=eval(['oo',int2str(i),'_.',jrc,'irfs.',char(var_list{j+(ifig-1)*9}),'_',char(strcat(shockid,varargin{2}(i))),';']);
                    else
                        IRF=eval(['oo',int2str(i),'_.',jrc,'irfs.',char(strcat(var_list{j+(ifig-1)*9},varargin{2}(i))),'_',char(strcat(shockid,varargin{2}(i))),';']);
                    end
                    
                end
            hold on
            plot(1:periods,flip*IRF(1:periods),'Linewidth',1.5,'MarkerSize',6)
            title(var_list{j+(ifig-1)*9},'interpreter','none')
            end
            plot(1:periods,zeros(periods,1),'Linewidth',1.5,'Color', [.7 .7 .7])
        end
        legend(lgnd,'interpreter','none')
    end
    if strcmp(color,'BW')
        set(groot,'defaultAxesColorOrder','remove')
        set(groot,'defaultAxesLineStyleOrder','remove')
    end
end


if strmatch('post_dens', comp_type, 'exact')
    if nargin==4
        nfig=ceil(length(var_list)/9);
            for ifig=1:nfig
                if ifig<nfig
                    nplots=9;
                else
                    nplots= rem(length(var_list),9);
                end
                figure(ifig)
                
                for j=1:nplots
                    subplot(3,3,j)
                    for i=1:length(folders)
                        
                        if strncmp(var_list{j+(ifig-1)*9},'EPS',3)==1
                            
                            
                            eval(['plot(oo',int2str(i),'_.posterior_density.shocks_std.',char(strcat(var_list{j+(ifig-1)*9},varargin{1}(i))),'(:,1),oo',int2str(i),'_.posterior_density.shocks_std.',char(strcat(var_list{j+(ifig-1)*9},varargin{1}(i))),'(:,2),''Linewidth'',2);']);
                            hold on
                        else
                            
                            
                            eval(['plot(oo',int2str(i),'_.posterior_density.parameters.',char(strcat(var_list{j+(ifig-1)*9},varargin{1}(i))),'(:,1),oo',int2str(i),'_.posterior_density.parameters.',char(strcat(var_list{j+(ifig-1)*9},varargin{1}(i))),'(:,2),''Linewidth'',2);']);
                            hold on
                        end
                        
                    end
                    if strncmp(var_list{j+(ifig-1)*9},'EPS',3)==1
                        
                        
                        eval(['plot(oo',int2str(i),'_.prior_density.shocks_std.',char(strcat(var_list{j+(ifig-1)*9},varargin{1}(i))),'(:,1),oo',int2str(i),'_.prior_density.shocks_std.',char(strcat(var_list{j+(ifig-1)*9},varargin{1}(i))),'(:,2),''Linewidth'',2,''Color'', [0.4 0.4 0.4]);']);
                        hold on
                    else
                        eval(['plot(oo',int2str(i),'_.prior_density.parameters.',char(strcat(var_list{j+(ifig-1)*9},varargin{1}(i))),'(:,1),oo',int2str(i),'_.prior_density.parameters.',char(strcat(var_list{j+(ifig-1)*9},varargin{1}(i))),'(:,2),''Linewidth'',2,''Color'', [0.4 0.4 0.4]);']);
                        hold on
                    end
                    title(var_list{j+(ifig-1)*9},'interpreter','none')
                end
                legend(lgnd,'interpreter','none')
            end
    else
        
        if ~ isequal(var_list ,'all')
            nfig=ceil(length(var_list)/9);
            for ifig=1:nfig
                if ifig<nfig
                    nplots=9;
                else
                    nplots= rem(length(var_list),9);
                end
                figure(ifig)
                
                for j=1:nplots
                    subplot(3,3,j)
                    for i=1:length(folders)
                        
                        if strncmp(var_list{j+(ifig-1)*9},'EPS',3)==1
                            
                            
                            eval(['plot(oo',int2str(i),'_.posterior_density.shocks_std.',var_list{j+(ifig-1)*9},'(:,1),oo',int2str(i),'_.posterior_density.shocks_std.',var_list{j+(ifig-1)*9},'(:,2),''Linewidth'',2);']);
                            hold on
                        else
                            eval(['plot(oo',int2str(i),'_.posterior_density.parameters.',var_list{j+(ifig-1)*9},'(:,1),oo',int2str(i),'_.posterior_density.parameters.',var_list{j+(ifig-1)*9},'(:,2),''Linewidth'',2);']);
                            hold on
                        end
                        
                    end
                    if strncmp(var_list{j+(ifig-1)*9},'EPS',3)==1
                        
                        
                        eval(['plot(oo',int2str(i),'_.prior_density.shocks_std.',var_list{j+(ifig-1)*9},'(:,1),oo',int2str(i),'_.prior_density.shocks_std.',var_list{j+(ifig-1)*9},'(:,2),''Linewidth'',2,''Color'', [0.4 0.4 0.4]);']);
                        hold on
                    else
                        eval(['plot(oo',int2str(i),'_.prior_density.parameters.',var_list{j+(ifig-1)*9},'(:,1),oo',int2str(i),'_.prior_density.parameters.',var_list{j+(ifig-1)*9},'(:,2),''Linewidth'',2,''Color'', [0.4 0.4 0.4]);']);
                        hold on
                    end
                    title(var_list{j+(ifig-1)*9},'interpreter','none')
                end
                legend(lgnd,'interpreter','none')
            end
        else
            clear var_list
            %parameters
            var_list=[];
            for i=1:length(folders)
                paramsize(i)=eval(['length(fieldnames(oo',int2str(i),'_.posterior_density.parameters));']);
                var_list=[var_list; eval(['fieldnames(oo',int2str(i),'_.posterior_density.parameters);'])];
            end
            var_list=unique(var_list);
            [nfig imaxfold]=max(paramsize);
            nfig=ceil(nfig/9);
            
            for ifig=1:nfig
                if ifig<nfig
                    nplots=9;
                else
                    nplots= rem(length(var_list),9);
                end
                figure(ifig)
                for j=1:nplots
                    subplot(3,3,j)
                    for i=1:length(folders)
                        
                        
                        
                        if eval(['isfield(oo',int2str(i),'_.posterior_density.parameters,''',var_list{j+(ifig-1)*9},''')'])
                            eval(['p(',int2str(i),')=plot(oo',int2str(i),'_.posterior_density.parameters.',var_list{j+(ifig-1)*9},'(:,1),oo',int2str(i),'_.posterior_density.parameters.',var_list{j+(ifig-1)*9},'(:,2),''Linewidth'',2);']);
                            hold on
                            ax = gca;
                            ax.ColorOrderIndex = ax.ColorOrderIndex-1;
                            eval(['plot(oo',int2str(i),'_.prior_density.parameters.',var_list{j+(ifig-1)*9},'(:,1),oo',int2str(i),'_.prior_density.parameters.',var_list{j+(ifig-1)*9},'(:,2),'':'',''Linewidth'',1);']);
                            hold on
                        else
                            ax = gca;
                            ax.ColorOrderIndex = ax.ColorOrderIndex+1;
                        end
                        hold on
                        
                        
                    end
                    %                 if strncmp(var_list{j+(ifig-1)*9},'EPS',3)==1
                    %
                    %
                    %                     eval(['plot(oo',int2str(imaxfold),'_.prior_density.shocks_std.',var_list{j+(ifig-1)*9},'(:,1),oo',int2str(imaxfold),'_.prior_density.shocks_std.',var_list{j+(ifig-1)*9},'(:,2),''Linewidth'',2,''Color'', [0.4 0.4 0.4]);']);
                    %                     hold on
                    %                 else
                    %                     eval(['plot(oo',int2str(imaxfold),'_.prior_density.parameters.',var_list{j+(ifig-1)*9},'(:,1),oo',int2str(imaxfold),'_.prior_density.parameters.',var_list{j+(ifig-1)*9},'(:,2),''Linewidth'',2,''Color'', [0.4 0.4 0.4]);']);
                    %                     hold on
                    %                 end
                    title(var_list{j+(ifig-1)*9},'interpreter','none')
                end
                legend(p,lgnd,'interpreter','none')
            end
            
            %shock_std
            ifigstart=ifig+1;
            for i=1:length(folders)
                shocksize(i)=eval(['length(fieldnames(oo',int2str(i),'_.posterior_density.shocks_std));']);
            end
            [nshock imaxfold]=max(shocksize);
            nfig=nfig+ceil(nshock/9);
            var_list=eval(['fieldnames(oo',int2str(imaxfold),'_.posterior_density.shocks_std);']);
            for ifig=ifigstart:nfig
                if ifig<nfig
                    nplots=9;
                else
                    nplots= rem(length(var_list),9);
                end
                figure(ifig)
                for j=1:nplots
                    subplot(3,3,j)
                    for i=1:length(folders)
                        
                        
                        if eval(['isfield(oo',int2str(i),'_.posterior_density.shocks_std,''',var_list{j+(ifig-ifigstart)*9},''')'])
                            
                            eval(['p(',int2str(i),')=plot(oo',int2str(i),'_.posterior_density.shocks_std.',var_list{j+(ifig-ifigstart)*9},'(:,1),oo',int2str(i),'_.posterior_density.shocks_std.',var_list{j+(ifig-ifigstart)*9},'(:,2),''Linewidth'',2);']);
                            hold on
                            ax = gca;
                            ax.ColorOrderIndex = ax.ColorOrderIndex-1;
                            eval(['plot(oo',int2str(i),'_.prior_density.shocks_std.',var_list{j+(ifig-ifigstart)*9},'(:,1),oo',int2str(i),'_.prior_density.shocks_std.',var_list{j+(ifig-ifigstart)*9},'(:,2),'':'',''Linewidth'',1);']);
                        else
                            ax = gca;
                            ax.ColorOrderIndex = ax.ColorOrderIndex+1;
                        end
                        hold on
                        
                        
                    end
                    
                    
                    %                 if strncmp(var_list{j+(ifig-ifigstart)*9},'EPS',3)==1
                    %
                    %
                    %                     eval(['plot(oo',int2str(imaxfold),'_.prior_density.shocks_std.',var_list{j+(ifig-ifigstart)*9},'(:,1),oo',int2str(imaxfold),'_.prior_density.shocks_std.',var_list{j+(ifig-ifigstart)*9},'(:,2),''Linewidth'',2,''Color'', [0.4 0.4 0.4]);']);
                    %                     hold on
                    %                 else
                    %                     eval(['plot(oo',int2str(imaxfold),'_.prior_density.parameters.',var_list{j+(ifig-ifigstart)*9},'(:,1),oo',int2str(imaxfold),'_.prior_density.parameters.',var_list{j+(ifig-ifigstart)*9},'(:,2),''Linewidth'',2,''Color'', [0.4 0.4 0.4]);']);
                    %                     hold on
                    %                 end
                    title(var_list{j+(ifig-ifigstart)*9},'interpreter','none')
                end
                legend(p,lgnd,'interpreter','none')
            end
        end
    end
end

if   strmatch('shock_decomp', comp_type, 'exact')
    nplots=length(folders);
    indx = strfind(var_list,'_');  
%     ootemp_=oo1_;
%     optionstemp_=options1_;
%     Mtemp_=M1_;
    for i=1:nplots
       eval(['options',int2str(i),'_.plot_shock_decomp.interactive=0;']);
       eval(['options',int2str(i),'_.plot_shock_decomp.detail_plot = 0;']);
        eval(['options',int2str(i),'_.plot_shock_decomp.use_shock_groups = var_list(indx+1:end);']);
        eval(['options',int2str(i),'_.plot_shock_decomp.type = ''aoa'';']);
        eval(['options',int2str(i),'_.plot_shock_decomp.nodisplay = 0;']);
        eval(['options',int2str(i),'_.plot_shock_decomp.write_xls = 0;']);
        eval(['options',int2str(i),'_.plot_shock_decomp.realtime = 1;']);
        eval(['options',int2str(i),'_.plot_shock_decomp.plot_init_date = dates(''2009Q1'');']);
         eval(['options',int2str(i),'_.plot_shock_decomp.plot_end_date = dates(''2019Q4'');']);
        eval(['options',int2str(i),'_.plot_shock_decomp.vintage = 0;']);
        eval(['options',int2str(i),'_.nobs = 89;']);
        eval(['assignin(''base'',''options_'',options',int2str(i),'_)']);
       eval(['[z',int2str(i),' shock_names my_initial_date]=compare_plot_shock_decomp(M',int2str(i),'_,oo',int2str(i),'_,options',int2str(i),'_,var_list);']); 
    end
    j=0;
    for i=1:nplots^2
        ind=1:nplots+1:nplots^2;
        
        
        if any(ind==i)
           j=j+1; 
        subplot(nplots,nplots,i)
        title(lgnd{j},'interpreter','none')
        eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.detail_plot = 0;']);
        eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.use_shock_groups = var_list(indx+1:end);']);
       eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.type = ''aoa'';']);
        eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.nodisplay = 0;']);
        eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.write_xls = 0;']);
         eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.realtime = 1;']);
            eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.plot_init_date = dates(''2009Q1'');']);
         eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.plot_end_date = dates(''2019Q4'');']);
                eval(['options',int2str(i-nplots*(j-1)),'_.plot_shock_decomp.vintage = 0;']);
                 eval(['options',int2str(i-nplots*(j-1)),'_.nobs = 89;']);
                  eval(['assignin(''base'',''options_'',options',int2str(i-nplots*(j-1)),'_)']);
        eval(['compare_plot_shock_decomp(M',int2str(i-nplots*(j-1)),'_,oo',int2str(i-nplots*(j-1)),'_,options',int2str(i-nplots*(j-1)),'_,var_list);']);
       
        else
%             [i_var,nvar] = varlist_indices(var_list,M1_.endo_names);
            subplot(nplots,nplots,i)
            if i<=(nplots*j)
                title([int2str(i-nplots*(j-1)),'-',int2str(j)])
                eval(['z=z',int2str(i-nplots*(j-1)),'-z',int2str(j),';']);
                compare_graph_decomp(z,shock_names,M1_.endo_names,1,dates('2009Y'),M1_,options1_);
                
                % eval(['ootemp_.shock_decomposition(i_var,:,:) = oo',int2str(i-nplots*(j-1)),'_.shock_decomposition(i_var,:,:)-oo',int2str(j),'_.shock_decomposition(i_var,:,:);']);
            else
                title([int2str(i-nplots*j),'-',int2str(j+1)])
                eval(['z=z',int2str(i-nplots*j),'-z',int2str(j+1),';']);
                compare_graph_decomp(z,shock_names,M1_.endo_names,1,dates('2009Y'),M1_,options1_);
                
                %eval(['ootemp_.shock_decomposition(i_var,:,:) = oo',int2str(i-nplots*j),'_.shock_decomposition(i_var,:,:)-oo',int2str(j+1),'_.shock_decomposition(i_var,:,:);']);
            end
       
       
        end
    end
end


if strmatch('data', comp_type, 'exact')
  for i=1:length(folders)
                var_list=[];
                var_list=[var_list; eval(['options',int2str(i),'_.varobs;'])];
               % eval(['load dataobs',int2str(i),'']) 

    for j=1:size(var_list,2)
    eval([char(strcat(var_list{1,j},'_',int2str(i))),'=dataobs',int2str(i),'.',var_list{1,j};]);
    eval(['clear ' char(var_list{j})])
    end
  end
  for j=1:size(var_list,2)
    figure(j),
    for i=1:length(folders)
            eval(['plot(',char(strcat(var_list{1,j},'_',int2str(i))),',''LineWidth'',2)']);
            hold on
    end
    hold off
    title(var_list{j},'interpreter','none')
    legend(lgnd,'interpreter','none')

  end
      
end

if  strmatch('mode', comp_type, 'exact')
    if nargin==4
        sizemode=varargin{1};
    else
        sizemode=20;
    end
    xparam_diff = abs( (xparam1_1 - xparam1_2)./xparam1_2 );
    [xparam_diff_sort, idx] = sort(xparam_diff,'descend');
    parameter_names_sort=parameter_names_2(idx);
    xparam1_1 = xparam1_1(idx);
    xparam1_2 = xparam1_2(idx);
    c=categorical(parameter_names_sort(1:sizemode));
    bar(c,xparam_diff_sort(1:sizemode))
    
    table(parameter_names_sort(1:sizemode),xparam1_1(1:sizemode),xparam1_2(1:sizemode))
    
end

