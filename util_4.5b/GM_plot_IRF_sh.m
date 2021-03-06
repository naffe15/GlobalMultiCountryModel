
function GM_plot_IRF_sh(xname, vargin, s1, texname, nr, nc, dirname, labels, Tend, init)
global M_ oo_ options_

fname_ =M_.fname;

if nargin<4 , texname=''; end
if nargin<7 || isempty(dirname), dirname='Output'; end
if nargin<8 || isempty(labels), xscale=1; end
if ~isfield(M_,'dname'),
    M_.dname = M_.fname;
end
DirectoryName = CheckPath(dirname,M_.dname);
% DirectoryName =pwd;
if nargin<5,
    nr = 3;
    nc = 4;
end

nsub = nc*nr;
if nargin<9 || isempty(Tend),
    Tend = options_.irf;
else
    Tend = min(options_.irf, Tend);
end

if nargin<10 || isempty(init),
    init = 1;
end

if init,
    if ~ismac
        delete([DirectoryName '\GM_IRF.TeX']);
    else
        delete([DirectoryName '/GM_IRF.TeX']);
    end
end


if nargin==0,
    return
end
if ~iscell(vargin)
    for j=1:size(vargin,1);
        vargin0{j,1}=deblank(vargin(j,:));
    end
    vargin=vargin0;
end
texnameflag =0;
if size(vargin,2)==3,
    texname=vargin(:,3);
    vargin=vargin(:,1:2);
    texnameflag =1;
end
if ~isempty(texname),
    texnameflag =1;
end    
if size(vargin,2)==2,
    temp1=vargin(:,2);
    vargin=vargin(:,1);
    for j=1:length(vargin),
        vscale(j,1)=temp1{j,1};
    end
else
    vscale=ones(length(vargin),1);
end
% Names of the shocks
if ~iscell(xname)
    for j=1:size(xname,1);
        xname0{j}=deblank(xname(j,:));
    end
    xname=xname0;
end
%
iexo=[];
if size(xname,2)==2,
    exoscale=cell2mat(xname(:,2));
    xname=xname(:,1);
else
    exoscale=ones(size(xname,1),1);
end
for i=1:size(xname,1)
    itemp = strmatch(xname{i},M_.exo_names,'exact');
    if isempty(itemp)
        error(['Shock ',xname{i},' is not defined!'])
    else
        iexo=[iexo, itemp];
    end
end
xname=cellstr(M_.exo_names(iexo,:));
% List of countries 
iendo=[];countries={};
ic=1;
for i=1:length(vargin),
    itemp = strmatch(vargin{i},M_.endo_names,'exact');
    if isempty(itemp)
        error(['Variable ',vargin{i},' is not defined!'])
    else        
        iendo=[iendo, itemp];
        k=strfind(vargin{i},'_');
        if ~isempty(k)
            countries=[countries;vargin{i}(k(end)+1:end)];
            ic=ic+1;
        end
    end
end
countries=unique(countries);
% For coloring place RoW for the first item
i=strmatch('RoW',countries);
if isempty(i)
    % Add RoW is not exist 
    countries=['RoW';countries];
else
    countries=[countries(i,:);countries(1:i-1,:);countries(i+1:end,:)];
end
vargin=cellstr(M_.endo_names(iendo,:));
% Find variable names without countrycode
varswithoutco=[];
isglobalvariable=0;
coindxs=zeros(length(vargin),length(countries));
for i=1:length(vargin)
    k=strfind(vargin{i},'_');
    if isempty(k) % Global variable
        varswithoutco{i}=vargin{i};
        isglobalvariable=1;
    else
        varswithoutco{i}=vargin{i}(1:k(end)-1);
        j=1;
        tempi=strfind(vargin{i}(k(end)+1:end),countries{j});
        while isempty(tempi)
            j=j+1;
            tempi=strfind(vargin{i}(k(end)+1:end),countries{j});
        end
        coindxs(i,j) = 1;
    end
end
% Find unique ones and keep the same order as in vargin
[uvars,ia,ic]=unique(varswithoutco,'stable');
vscale=vscale(ia);
% For each uvars find the corresponding countries
uvarsco=zeros(length(uvars),length(countries));
for i=1:length(uvars)
    j=find(ic==i);
    uvarsco(i,:)=sum(coindxs(j,:),1);
end
uvarsco(find(uvarsco>1))=1;
% nv=length(vargin);
nv=length(uvarsco);
%
%
if (nargin<4 && options_.TeX) || (isempty(texname) && options_.TeX)
    for j=1:length(vargin),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end
func = @(x) colorspace('RGB->Lab',x);
% First color is for the global variable and the second for RoW
colors = distinguishable_colors(length(countries)+1,'w',func);  
c=colors(1,:);
colors(1,:)=colors(2,:);
colors(2,:)=c;

style={'-','--','-.'};
for i=1:length(iexo)  % For all shocks
    nfig=0;
    nplo=0;
    indplot=[];
    for j=1:length(uvars) % For all variables
        globalvariable=0;
        oktoplot=0;
        if sum(uvarsco(j,:))>0
            % Variable with countrycode
            for c=1:size(uvarsco,2) % For all countries
                if uvarsco(j,c)==1,
                    name = [uvars{j} '_' countries{c} '_' deblank(M_.exo_names(iexo(i),:))];
                    if texnameflag,
                        inx = strmatch([uvars{j} '_' countries{c}],vargin);
                        tname = texname{inx};
                    end
                    try
                        eval(['MeanIRF=s1.' name,';']);
                        %eval(['MeanIRF=rand(1,Tend);']);
                    catch
                        MeanIRF = 0;
                    end
                    if max(abs(MeanIRF)) > 1e-6
                        oktoplot=1;
                    end
                end
            end
        else
            % Global variable
            globalvariable=1;
            name = [uvars{j} '_' deblank(M_.exo_names(iexo(i),:))];
            if texnameflag,
                inx = strmatch(uvars{j},vargin);
                tname = texname{inx};
            end
            try
                eval(['MeanIRF=s1.' name,';']);
                %eval(['MeanIRF=rand(1,Tend);']);
            catch
                MeanIRF = 0;
            end
            if max(abs(MeanIRF)) > 1e-6
                oktoplot=1;
            end
        end
        if oktoplot==1
            nplo=nplo+1;
            indplot=[indplot;j];
            if globalvariable
                fname=uvars{j};                
            else
                fname=[uvars{j} '_' '@{co}'];        
            end
            if texnameflag,
                fname = tname;
            end
            if mod(nplo,nsub)==1
                dyn_figure(options_.nodisplay,'name',['Orthogonalised shocks to ',deblank(M_.exo_names(iexo(i),:))]);
                % figure('name',['Orthogonalised shocks to ',deblank(M_.exo_names(iexo(i),:))]);
                nfig=nfig+1;
                nplo = 1;
            end
            subplot(nr,nc,nplo);
            if nargin==8,
                xscale=labels{j,3};
            end
            plot([1 Tend]/xscale,[0 0],':k','linewidth',0.5);
            hold on;
            if globalvariable==1
                name = [uvars{j} '_' deblank(M_.exo_names(iexo(i),:))];
                try
                    eval(['MeanIRF=s1.' name,';']);
                    %eval(['MeanIRF=(-0.5:1/(Tend-1):0.5);']);
                catch
                    MeanIRF = 0;
                end
                if max(abs(MeanIRF)) > 1e-6
                    plot([1:Tend]/xscale,MeanIRF(1:Tend)*vscale(j)*exoscale(i),'color',colors(1,:),'linewidth',1);
                end
            else
                for c=1:size(uvarsco,2) % For all countries
                    if uvarsco(j,c)==1
                        name = [uvars{j} '_' countries{c} '_' deblank(M_.exo_names(iexo(i),:))];
                        try
                            eval(['MeanIRF=s1.' name,';']);
                            % eval(['MeanIRF=c*ones(1,Tend);']);
                        catch
                            MeanIRF = 0;
                        end
                        if max(abs(MeanIRF)) > 1e-6
                            plot([1:Tend]/xscale,MeanIRF(1:Tend)*vscale(j)*exoscale(i),'color',colors(c+1,:),'linewidth',1.5,'Linestyle',style{c});
                        end
                    end
                end
            end
            xlim([1 Tend]/xscale);
            if nargin==8,
                text(0.82,-0.2,labels{j,1},'units','normalized');
                text(0.03,0.89,labels{j,2},'units','normalized');
            end
            if options_.TeX==0,                
                %title(texname{j},'interpreter','tex')
                title(fname,'interpreter','none')
            else
%                 title(tname,'interpreter','none')
                title(fname,'interpreter','none')
            end
        end
            if (mod(nplo,nsub)==0 || j==length(uvars)) && nplo,
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
                for k = 1:p
                    h = axes('position',[fL(k),0,ffl,0.05]);
                    %annotation('textbox', [fL(k),0,ffl,0.05],'String', country_code{k},'Color',colorsA(k,:),'horizontalalignment','center');
					annotation('line',[fL(k)+0.15 fL(k)+0.2 ],[0.025 0.025 ],'LineStyle',style{k},'Color',colorsA(k,:),'linewidth',1.5);
                    annotation('textbox', [fL(k),0,ffl,0.05],'String', country_code{k} ,'Color',colorsA(k,:),'horizontalalignment','center','verticalalignment','middle');

                    set(h,'YTick',[]);set(h,'XTick',[]);
                end               
                % End Annotation of the plot 
                if ~ismac
                dyn_saveas(gcf,[DirectoryName '\GM_IRF_',int2str(Tend),'_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
                else
                dyn_saveas(gcf,[DirectoryName '/GM_IRF_',int2str(Tend),'_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
                end
                if options_.TeX,
                    texfilename=[DirectoryName '\GM_IRF.TeX'];
                    l1=(nc*nr)*(nfig-1)+1;
                    
                    %         texname0=texname(posi,:);
                    %         fidTeX = fopen([M_.fname '_Fit_Plots.TeX'],'w+');
                    %         files=dir('*_Fit*.eps');
                    fidTeX = fopen([DirectoryName '\GM_IRF.TeX'],'a+');
                    %         tempfilename=[DirectoryName '\' fname_,'_IRF_',deblank(M_.exo_names(iexo(i),:)),'_1'];
                    %         tempfilename=deblank(tempfilename);
                    %         [a,b,c]=fileparts(tempfilename);
                    %         ab=strcat(a,'/',b);
                    ab=([DirectoryName '/GM_IRF_',int2str(Tend),'_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
                    fprintf(fidTeX,'\\begin{figure}[H]\n');
                    %for kk=1:nplo,
                    %    fprintf(fidTeX,'\\psfrag{%s}[1][][1][0]{%s}\n',deblank(vargin{indplot(l1)}),['$', deblank(texname{indplot(l1)}), '$']);
                    %    l1=l1+1;
                    %end
                    %         numberofvars=numberofvars-9;
                    fprintf(fidTeX,'\\centering \n');
                    fprintf(fidTeX,['\\includegraphics[scale=0.7] {%s} \n'],ab);
                    fprintf(fidTeX,'\\caption{Impulse responses to a shock to $%s$}',[deblank(M_.exo_names_tex(iexo(i),:))]);
                    fprintf(fidTeX,'\\label{Fig: IRF Plot:%s:%s}\n',int2str(Tend),int2str(j));
                    fprintf(fidTeX,'\\end{figure}\n');
                    fprintf(fidTeX,' \n');
                    fclose(fidTeX);
                end
                nplo=0;
            end
    end
end

%         movefile(texfilename,DirectoryName);
