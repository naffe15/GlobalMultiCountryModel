function GM_plot_smooth(varargin)
global options_ M_ oo_
% persistent h;

ys_ = [oo_.steady_state*0; zeros(M_.exo_nbr,1)];
lgy_ = char(M_.endo_names,M_.exo_names);
fname_ = M_.fname;
texname=char(M_.endo_names_tex,M_.exo_names_tex);
SmoothedVariables=[struct2cell(oo_.SmoothedVariables); struct2cell(oo_.SmoothedShocks)];
my_field_names = [fieldnames(oo_.SmoothedVariables); fieldnames(oo_.SmoothedShocks)];
isvar=zeros(length(SmoothedVariables),1);
for jf = 1:length(SmoothedVariables),
    isvar(jf)=~(isstruct(SmoothedVariables{jf}));
end
SmoothedVariables=cell2struct(SmoothedVariables(logical(isvar)),my_field_names(logical(isvar)));

ifig0=0;ivar=[];countries={};
cooption=0;vararginorig={};
if strmatch('-',varargin{end})
    while strmatch('-',varargin{end})
        optn = varargin{end};
        varargin=varargin(1:end-1);
        if strmatch('-new',optn)
            ifig0=length(dir(['GM_Smooth*.fig']));
        end
        if strmatch('-co',optn)            
            cooption=1;
            % There is country option 
            i1=strfind(optn,'(');
            i2=strfind(optn,')');
            % syntax -co() or -co means that use all countries 
            if (isempty(i1) && isempty(i2)) || (i2-i1)==1
                % Use all possible countries
                svnames=fieldnames(SmoothedVariables);
                for i=1:length(varargin)
                    if isempty(strmatch('-new',varargin{i}))                        
                        k1=strmatch(varargin{i},svnames);
                        kco=[];kgo=[];
                        for k=1:length(k1)
                            s=char(svnames(k1(k)));
                            if length(strfind(s(length(varargin{i})+1:end),'_'))==1 & strfind(s(length(varargin{i})+1),'_')                                                                    
                                kco=[kco;k1(k)]; % Variable with a country code
                            elseif length(s)==length(varargin{i})                                
                                kgo=[kgo;k1(k)];  % Global variables
                            end
                        end                                      
                        if isempty(kco) & isempty(kgo)
                            error(['Variable ',varargin{i},' is not defined!'])
                        end                     
                        if isempty(kco)
                            % Global variable
                            itemp = strmatch(varargin{i},fieldnames(SmoothedVariables),'exact');
                            ivar=[ivar, itemp];
                        else
                            k=kco;
                            % svnames(k) contains the full names
                            for j=1:length(k),
                                itemp = strmatch(svnames{k(j)},fieldnames(SmoothedVariables),'exact');
                                ivar=[ivar, itemp];
                                k1=strfind(svnames{k(j)},'_');
                                if ~isempty(k1)
                                    countries=[countries;svnames{k(j)}(k1(end)+1:end)];
                                end
                            end
                        end
                    end
                    %if ~isempty(k)
                    %    countries=[countries;varargin{i}(k(end)+1:end)];
                    %    ic=ic+1;
                    %end
                end
            else
                % Use countries X,Y,Z,... separated by commas in -co(X,Y,Z,...)
                co=[];
                k=i1+1;
                while k<=(i2-1)
                    if optn(k)~=','
                        co=[co optn(k)];
                    end
                    if optn(k)==',' || k==(i2-1)
                        countries=[countries;co];
                        co=[];
                    end
                    k=k+1;
                end
            end
            countries=unique(countries);
            % Find full names of the variables
            svnames=fieldnames(SmoothedVariables);         
            cofound=zeros(length(countries),1);
            for i=1:length(varargin)
                vname=varargin{i};
                if isempty(strmatch('-new',vname))
                    ii=[strfind(vname,'[') strfind(vname,'(')];
                    if ~isempty(ii)
                        vname(ii:end)=[];
                    end
                    k1=strmatch(vname,svnames);
                    k2=strmatch([vname '_'],svnames);
                    if isempty(k1) && isempty(k2)
                        error(['Variable ',vname,' is not defined!'])
                    else
                        k=[];
                        if ~isempty(k2)
                            if length(k1)~=length(k2)
                                k=k2;
                            else
                                k=k1;
                            end
                        end
                    end
                    if isempty(k)
                        % Global variable
                        itemp = strmatch(vname,fieldnames(SmoothedVariables),'exact');
                        ivar=[ivar itemp];
                        if ~isempty(ii)
                            vararginorig=[vararginorig;[vname varargin{i}(ii:end)]];
                        else
                            vararginorig=[vararginorig;vname];
                        end
                    else
                        for c=1:length(countries)
                            itemp=strmatch([vname '_' countries{c}],fieldnames(SmoothedVariables),'exact');
                            if ~isempty(itemp)
                                cofound(c)=1;
                                ivar=[ivar itemp];
                                if ~isempty(ii)
                                    vararginorig=[vararginorig;[vname '_' countries{c} varargin{i}(ii:end)]];
                                else
                                    vararginorig=[vararginorig;[vname '_' countries{c}]];
                                end
                            end
                        end
                    end
                end
            end
            if matlab_ver_less_than('8')
                ivar=unique(ivar);
            else
                ivar=unique(ivar,'stable');
            end
            % Produce an error message if there is country that
            % is not found at all from listed variables
            c=find(cofound==0);
            if ~isempty(c)
                notdefc=[];
                for k=1:length(c)
                    notdefc=[notdefc countries{c(k)} ' '];
                end
                error(['Country/countries ',notdefc,' not defined!'])
            end
        end
    end
    v=fieldnames(SmoothedVariables);
    varargin=cellstr(v(ivar,:))';
end
% For coloring place RoW for the first item
i=strmatch('RoW',countries);
if isempty(i)
    % Add RoW is not exist 
    countries=['RoW';countries];
else
    countries=[countries(i,:);countries(1:i-1,:);countries(i+1:end,:)];
end
if ifig0==0
    afig=dir(['GM_Smooth*.*']);
    for j=1:length(afig),
        delete(afig(j).name);
    end
    delete(['GM_Smooth_Plots.TeX'])
end
if exist(options_.datafile)
    instr = options_.datafile;
else
    instr = ['load ' options_.datafile];
end
try
    eval(instr);
    if ~exist('T','var'),
        temp = eval(deblank(options_.varobs(1,:)));
        T=[1:length(temp)]';
        clear temp;
    end
catch
    T=[1:options_.nobs+options_.first_obs-1];
end

fobs = options_.first_obs;
if cooption==0
    % Find list of countries
    countries={};ic=1;
    for i=1:length(varargin),
        k=strfind(varargin{i},'_');
        if ~isempty(k)
            countries=[countries;varargin{i}(k(end)+1:end)];
            ic=ic+1;
        end
    end
    countries=unique(countries);
    vararginorig=varargin;
end
% Find variable names without countrycode
varswithoutco=[];
isglobalvariable=0;
coindxs=zeros(length(varargin),length(countries));
for i=1:length(varargin)    
    k=strfind(varargin{i},'_');
    if isempty(k) % Global variable 
        varswithoutco{i}=varargin{i};
        isglobalvariable=1;
    else
        varswithoutco{i}=varargin{i}(1:k(end)-1);
        j=1;
        tempi=strcmp(varargin{i}(k(end)+1:end),countries{j});
        while ~tempi
            j=j+1;
            tempi=strcmp(varargin{i}(k(end)+1:end),countries{j});
        end
        coindxs(i,j) = 1;        
    end
end
% Find unique ones and keep the same order as in vargin
if matlab_ver_less_than('8')
    [uvars,ia,ic]=unique(varswithoutco);
else
    [uvars,ia,ic]=unique(varswithoutco,'stable');
end
% vscale=vscale(ia);
% For each uvars find the corresponding countries
for i=1:length(uvars)
    j=find(ic==i);
    uvarsco(i,:)=sum(coindxs(j,:),1);
end
uvarsco(find(uvarsco>1))=1;
if isempty(countries)
    uvarsco=zeros(length(uvars),1);
end

% nv=length(varargin);
nv=length(uvarsco);
func = @(x) colorspace('RGB->Lab',x);
% First color is for the global variable and the second for RoW
colors = distinguishable_colors(length(countries)+1,'w',func);  
c=colors(1,:);
colors(1,:)=colors(2,:);
colors(2,:)=c;
ttrend = [1:length(T(fobs:fobs+options_.nobs-1))]';
nfig=ifig0;
nplo=0;
nsub=9; 
for ii=1:length(uvars)    
    nplo=nplo+1;
    if mod(nplo,nsub)==1        
        h = dyn_figure(options_.nodisplay,'Name','Smooth variables');    
        nfig=nfig+1;        
        %figure(nfig);
        nplo=1;         
    end
    subplot(3,3,nplo);hold on;
    % All countries with the variable uvarsco{j}
    k=1;
    while k<=size(uvarsco,2)
        if sum(uvarsco(ii,:))==0
            globalvariable=1;
            name=[uvars{ii}];
            k=size(uvarsco,2)+1;
            colorindex=1; 
        else
            globalvariable=0;
            name=[uvars{ii} '_' countries{k}];
            colorindex=k+1;
            k=k+1;            
        end
        if uvarsco(ii,k-1) || globalvariable==1
            j=strmatch(name,varargin,'exact');
            exponential_trend=0;
            trendvar{j} = '';
            if ~isempty(strfind(vararginorig{j},'('))
                i0 = strfind(vararginorig{j},'(');
                i1 = strfind(vararginorig{j},')')-1;
                trend(j) = get_param_by_name(vararginorig{j}(i0+1:i1));
                %vararginorig{j}=vararginorig{j}(1:i0-1);
            elseif ~isempty(strfind(vararginorig{j},'['))
                i0 = strfind(vararginorig{j},'[');
                i1 = strfind(vararginorig{j},',')-1;
                if isempty(i1),
                    i1 = strfind(vararginorig{j},']')-1;
                else
                    i2 = strfind(vararginorig{j},']')-1;
                    % Check if trendvar is global or country specific
                    v1=strmatch(vararginorig{j}(i1+2:i2),fieldnames(SmoothedVariables),'exact');  % Global
                    v2=[];  % country specific
                    if uvarsco(ii,k-1)
                        v2=strmatch([vararginorig{j}(i1+2:i2) '_' countries{k-1}] ,fieldnames(SmoothedVariables),'exact');
                    end
                    if ~isempty(v1)
                        trendvar{j} = vararginorig{j}(i1+2:i2);                            
                    elseif ~isempty(v2)
                        trendvar{j} = [vararginorig{j}(i1+2:i2) '_' countries{k-1}];                         
                    else
                       error(['Trend variable ',vararginorig{j}(i1+2:i2),' or ',[vararginorig{j}(i1+2:i2) '_' countries{k-1}],' is not a defined!'])
                    end
                end
                %trend(j) = evalin('base',varargin{j}(i0+1:i1));
                % Check is the vararginorig{j}(i0+1:i1) is global or
                % country specific 
                v1=strmatch(vararginorig{j}(i0+1:i1),M_.param_names,'exact');  % Global 
                v2=[];  % country specific 
                if uvarsco(ii,k-1)
                    v2=strmatch([vararginorig{j}(i0+1:i1) '_' countries{k-1}] ,M_.param_names,'exact');
                end
                if ~isempty(v1) 
                    trend(j) = get_param_by_name(vararginorig{j}(i0+1:i1));
                elseif ~isempty(v2)
                    trend(j) = get_param_by_name([vararginorig{j}(i0+1:i1) '_' countries{k-1}]);
                else
                    error(['Trend parameter',vararginorig{j}(i0+1:i1),' or ',[vararginorig{j}(i0+1:i1) '_' countries{k-1}],' is not a defined!'])
                end
                %vararginorig{j}=vararginorig{j}(1:i0-1);
                exponential_trend=1;
            else
                trend(j) = 0;
            end
            i = strmatch(varargin{j},lgy_,'exact');
            posi(j)=i;
            if ~isempty(texname)
                if strmatch('\log',texname(i,:));
                    texname(i,1:end-1)=texname(i,2:end);
                end
            end
            if exist(varargin{j},'var')
                yobs = eval([varargin{j},'(fobs:fobs+options_.nobs-1)']);
            end
            if exponential_trend==1,
                if ~isempty(trendvar{j}),
                    iplus = strfind(trendvar{j},'+');
                    iminus = strfind(trendvar{j},'-');
                    if ~isempty(iplus),
                        techno=get_smooth(trendvar{j}(1:iplus(1)-1)) + get_smooth(trendvar{j}(iplus(1)+1:end));
                    elseif ~isempty(iminus),
                        techno=get_smooth(trendvar{j}(1:iminus(1)-1)) - get_smooth(trendvar{j}(iminus(1)+1:end));
                    else
                        techno = get_smooth(trendvar{j});
                    end
                else
                    techno=0;
                end
                ymodel = eval(['[(SmoothedVariables.',varargin{j},'(1:end)+ys_(i)).*exp(trend(j)*ttrend+techno)]']);
                if exist(varargin{j},'var') && trend(j)
                    ymodel = ymodel.*mean_nan(yobs./ymodel);
                end
            else
                ymodel = eval(['[SmoothedVariables.',varargin{j},'(1:end)+ys_(i)+trend(j)*ttrend]']);
                if exist(varargin{j},'var') && trend(j)
                    ymodel = ymodel+mean_nan(yobs-ymodel);
                end
            end
            plot(T(fobs:fobs+options_.nobs-1),ymodel,'Color',colors(colorindex,:))
            if exist(varargin{j},'var')
                hold on,
                %             eval(['plot(T(fobs:fobs+options_.nobs-1),',varargin{j},'(fobs:fobs+options_.nobs-1),''-r.'')'])
                plot(T(fobs:fobs+options_.nobs-1),yobs,':.','Color',colors(colorindex,:))
            end
            
            if options_.TeX
                %             set(get(gca,'title'),'string',varargin{j},'interpreter','none')
                % title(varargin{j},'interpreter','none')
                if globalvariable==1
                    title(uvars{ii},'interpreter','none');
                else
                    title([uvars{ii} '_@{co}'],'interpreter','none');
                end
            else
                %             set(get(gca,'title'),'string',deblank(texname(i,:)),'interpreter','tex')
                if globalvariable==1
                    % title(deblank(texname(i,:)),'interpreter','tex')
                    title(uvars{ii},'interpreter','none');
                else
                    % title(deblank(texname(i,:)),'interpreter','tex')
                    title([uvars{ii} '_@{co}'],'interpreter','none');
                end
            end
            set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
        end
        if ((nplo==nsub || ii==length(uvars)) && (k > size(uvarsco,2)))
            % Annotation of the plot
            colorsA=colors;
            if isglobalvariable==1                
                if sum(uvarsco(:,1))>0
                    country_code=['Global';countries];                    
                else
                    country_code=['Global';countries(2:end,:)];
                    colorsA(2,:)=[];
                end
            else
                if sum(uvarsco(:,1))>0
                    country_code=countries;
                    colorsA(1,:)=[];
                else
                    country_code=countries(2:end,:);
                    colorsA(1:2,:)=[];
                end                
            end
            p=length(country_code);
            bf = 0.1;
            ffs = 0.05/(p-1);
            ffl = (1-2*bf-0.05)/p;
            if p>1,
                fL = linspace(bf,1-bf+ffs,p+1);
            else
                fL = bf;
            end
            for kk = 1:p
                h = axes('position',[fL(kk),0,ffl,0.05]);
                annotation('textbox', [fL(kk),0,ffl,0.05],'String', country_code{kk},'Color',colorsA(kk,:),'horizontalalignment','center');
                set(h,'YTick',[]);set(h,'XTick',[]);
            end
            dyn_saveas(gcf,['GM_Smooth',int2str(nfig)],options_.nodisplay,options_.graph_format)
            if options_.TeX,
                h1=1;
                %numberofvars=nv;
                %texname0=texname(posi,:);
                fidTeX = fopen(['GM_Smooth_Plots.TeX'],'a+');
                tempfilename=['GM_Smooth',int2str(nfig)];
                [a,b,c]=fileparts(tempfilename);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                %for kk=1:min(numberofvars,9)
                %    fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(varargin{kk}),['$' deblank(texname0(h1,:)) '$']);
                %    h1=h1+1;
                %end
                %numberofvars=numberofvars-9;
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,['\\includegraphics[width=0.80\\textwidth] {' b '} \n']);
                fprintf(fidTeX,'\\caption{Smooth}');
                fprintf(fidTeX,'\\label{Fig: Smooth :%s}\n',int2str(nfig));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,' \n');
                fclose(fidTeX);
            end
        end
    end
end