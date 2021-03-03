function GM_plot_simul(vargin, texname, nr, nc, dirname, tlim, out, legtxt,simul_name,init)
global M_ options_ oo_
persistent simul_nbr

fname_ =M_.fname;

if nargin<10 || isempty(init),
    init = 1;
end
if nargin<5 || isempty(dirname), dirname='Output'; end
if ~isfield(M_,'dname'),
    M_.dname = M_.fname;
end
DirectoryName = CheckPath(dirname,M_.dname);
% DirectoryName =pwd;

if nargin<3,
    nr = 3;
    nc = 4;
end

nsub = nc*nr;


if ~iscell(vargin)
    for j=1:size(vargin,1);
        vargin0{j,1}=deblank(vargin(j,:));
    end
    vargin=vargin0;
end
if size(vargin,2)==2,
    flag_scale=1;
    temp1=vargin(:,2);
    vargin=vargin(:,1);
    for j=1:length(vargin),
        vscale(j,1)=temp1{j,1};
    end
else
    flag_scale=0;
    vscale=ones(length(vargin),1);
    
end
% List of countries
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
    % Add RoW if not exist
    countries=['RoW';countries];
else
    countries=[countries(i,:);countries(1:i-1,:);countries(i+1:end,:)];
end
vargin=cellstr(M_.endo_names(iendo,:));
%
% Find unique ones and keep the same order as in vargin
%if matlab_ver_less_than('8')
%    [uvars,ia,ic]=unique(varswithoutco);
%else
%    [uvars,ia,ic]=unique(varswithoutco,'stable');
%end
%
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
if matlab_ver_less_than('8')
    [uvars,ia,ic]=unique(varswithoutco);
else
    [uvars,ia,ic]=unique(varswithoutco,'stable');
end
% vscale=vscale(ia);
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
if (nargin<2) || (isempty(texname) )
    for j=1:length(vargin),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end
if init,
    if options_.TeX,
        delete([DirectoryName '\GM_simul.TeX']);
    end
    clear simul_nbr
    simul_nbr = 1;
else
    simul_nbr = simul_nbr + 1;
end
if nargin<9 || isempty(simul_name),
    simul_name = int2str(simul_nbr);
end


if nargin<6 || isempty(tlim),
    name = vargin{1};
    [z,zss]=dyn2vec(name);
    tlim=[1 length(z)];
end
nfig=0;
nplo=0;
indplot=[];
func = @(x) colorspace('RGB->Lab',x);
% First color is for the global variable and the second for RoW
colors = distinguishable_colors(length(countries)+1,'w',func);
c=colors(1,:);
colors(1,:)=colors(2,:);
colors(2,:)=c;
for j=1:length(uvars),  % One variable at time
    %
    nplo=nplo+1;
    indplot=[indplot;j];
    if mod(nplo,nsub)==1,
        dyn_figure(options_.nodisplay,'name',['Deterministic simulation: ', simul_name]);
        nfig=nfig+1;        
        nplo = 1;
    end
    subplot(nr,nc,nplo);hold on;
    c=1;
    while c<=size(uvarsco,2)
        if sum(uvarsco(j,:))==0
            % Global variable
            name = uvars{j};
            globalvariable=1;
            cindex=1;
            c=size(uvarsco,2)+1;            
        else
            % Variable with country code
            if uvarsco(j,c)==1
                name=[ uvars{j} '_' countries{c}];
            end
            globalvariable=0;
            cindex=c+1;
            c=c+1;
            
        end
        if uvarsco(j,c-1) || globalvariable==1
            [z,zss]=dyn2vec(name);
            if nargin >=7,
                yy=zeros(length(z),length(out));
                for jy=1:length(out),
                    yy(:,jy)=getfield(out(jy),name);
                end
            end
            if flag_scale,
                z=z-zss;
                zss=0;
            end
            plot(z*vscale(j),'color',colors(cindex,:),'linewidth',1)
            if ~isempty(yy),
                hold all, plot(yy*vscale(j),'linewidth',1),
                if ~isempty(legtxt) && nplo==1,
                    %             txt=strcat(legtxt{:});
                    %             txt='';
                    %             for jj=1:length(legtxt),
                    %                 txt=[txt,'(',int2str(jj),') = ',legtxt{jj},'; '];
                    %             end
                    %             text(1,1.4,txt,'units','normalized','HorizontalAlignment','center')
                    legend(legtxt,'orientation','horizontal','position',[0.25 0.02 0.5 0.03])
                end
            end
            hold on,
            plot(tlim,[zss zss],':k','linewidth',0.5);
            xlim(tlim);
            % hold off
            %
            if globalvariable
                fname=uvars{j};
            else
                fname=[uvars{j} '_' '@{co}'];
            end
        end
    end
    if options_.TeX==0,
        % title(texname{j},'interpreter','tex')
        title(fname,'interpreter','none')
    else
        title(fname,'interpreter','none')
    end
    if (mod(nplo,nsub)==0 || j==length(vargin)) && nplo,
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
            annotation('textbox', [fL(k),0,ffl,0.05],'String', country_code{k},'Color',colorsA(k,:),'horizontalalignment','center');
            set(h,'YTick',[]);set(h,'XTick',[]);
        end
        dyn_saveas(gcf,[DirectoryName '\GM_simul_',simul_name,'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
        if options_.TeX,
            texfilename=[DirectoryName '\GM_simul.TeX'];
            %l1=(nc+nr)*(nfig-1)+1;
            
            fidTeX = fopen([DirectoryName '\GM_simul.TeX'],'a+');
            ab=([DirectoryName '/GM_simul_',simul_name,'_',int2str(nfig)]);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            %for kk=1:nplo,
            %    fprintf(fidTeX,'\\psfrag{%s}[1][][1][0]{%s}\n',deblank(vargin{indplot(l1)}),['$', deblank(texname{indplot(l1)}), '$']);
            %    l1=l1+1;
            %end
            %         numberofvars=numberofvars-9;
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[scale=0.7] {%s} \n',ab);
            fprintf(fidTeX,'\\caption{Deterministic simulation %s:%s}', int2str(simul_nbr), int2str(nfig));
            fprintf(fidTeX,'\\label{Fig:DetSim:%s:%s}\n',simul_name,int2str(nfig));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
            fclose(fidTeX);
        end
        nplo=0;
    end
end

save([DirectoryName '\' fname_,'_simul_',simul_name],'M_','oo_')


%         movefile(texfilename,DirectoryName);
