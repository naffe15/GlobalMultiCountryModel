function AnnForecastPlot = annualized_plot_simul(vargin, q2a, texname, nr, nc, dirname, tlim, flag_scale, simul_name,init)
global M_ options_ 
persistent simul_nbr

fname_ =M_.fname;

if nargin<10 || isempty(init),
    init = 1;
end
if nargin<6 || isempty(dirname), dirname='Output'; end
if ~isfield(M_,'dname'),
    M_.dname = M_.fname;
end
DirectoryName = CheckPath(dirname,M_.dname);
% DirectoryName =pwd;

if nargin<4,
    nr = 3;
    nc = 4;
end

nsub = nc*nr;

if nargin<8 || isempty(flag_scale)
    flag_scale=0;
end
vscale=ones(length(vargin),1);

iendo=[];
for i=1:length(vargin)
    itemp = strmatch(vargin{i},M_.endo_names,'exact');
    if isempty(itemp),
        error(['variable ',vargin{i},' does not exist']);
    else        
        iendo=[iendo, strmatch(vargin{i},M_.endo_names,'exact')];
    end
end
if (nargin<3) || (isempty(texname) )
    for j=1:length(vargin),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end
if options_.TeX,
    texfilename=[DirectoryName filesep fname_ '_annual_simul.TeX'];
end
if init,
    if options_.TeX,
        delete(texfilename);
    end
    clear simul_nbr
    simul_nbr = 1;
else
    simul_nbr = simul_nbr + 1;    
end
if nargin<9 || isempty(simul_name),
    simul_name = int2str(simul_nbr);
end

numberofvars=length(vargin);
nfig=0;
nplo=0;
indplot=[];

if nargin<7 || isempty(tlim),
    name = vargin{1};
    [z,zss]=dyn2vec(name);
    tlim=[1 length(z)];
end    

for j=1:length(vargin),
    
    name = vargin{j};
    [z,zss]=dyn2vec(name);
    z=z(M_.maximum_exo_lag+1:tlim(2));
    
    zpreamble = repmat(zss,8,1);
    ztmp = [zpreamble; z]-zss; 

    if isstruct(q2a(j).aux)
        zaux_name = q2a(j).aux.y;
        [zaux,zauxss]=dyn2vec(zaux_name);
        zaux=zaux(M_.maximum_exo_lag+1:tlim(2));
        
        zpreamble = repmat(zauxss,8,1);
        ztmp_aux = [zpreamble; zaux]-zauxss;
        q2a(j).aux.y = ztmp_aux;
    end
    
    [za, zass, gza, gzass] = ...
        quarterly2annual(ztmp,zss,q2a(j).GYTREND0,q2a(j).type,q2a(j).islog,q2a(j).aux);
    
    if q2a(j).plot ==1
        ztmp=gza(2:end)+gzass;
        aname = q2a(j).gname;
        atexname = q2a(j).tex_gname;
        zssa = gzass;
    elseif q2a(j).plot == 2
        ztmp=za(2:end)+zass;
        aname = q2a(j).name;
        atexname = q2a(j).tex_name;
        zssa = zass;
    else
        error('choose level or growth rate')
    end
    
    za=ztmp;
    if flag_scale
        za=za-zssa;
        zssa=0;
    end
    nplo=nplo+1;
    indplot=[indplot;j];
    if mod(nplo,nsub)==1,
        dyn_figure(options_.nodisplay,'name',['Deterministic simulation: ', simul_name]);
        nfig=nfig+1;
    end
    subplot(nr,nc,nplo)
    plot(0:length(za)-1,za*vscale(j),'-k','linewidth',2)
    AnnForecastPlot(j).VarName = aname;
    AnnForecastPlot(j).TimeLineQA = [0:length(za)-1]';
    AnnForecastPlot(j).exogassmpt = za;

    hold on,
    plot([-0.5 length(za)-0.5],[zssa zssa],'r','linewidth',0.5);
    xlim([0 length(za)-1]); 
%     xtick = get(gca,'xtick');
%     xtick=xtick(find(mod(xtick,4)==0));
    xlabel('periods (yrs)')
    hold off
    if options_.TeX==0
        title(atexname,'interpreter','tex')
    else
        title(aname,'interpreter','none')
    end
    if (mod(nplo,nsub)==0 || j==length(vargin)) && nplo,
        dyn_saveas(gcf,[DirectoryName filesep fname_,'_annual_simul_',simul_name,'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
        if options_.TeX
            l1=(nc+nr)*(nfig-1)+1;
            
            fidTeX = fopen(texfilename,'a+');
            ab=([DirectoryName '/' fname_,'_annual_simul_',simul_name,'_',int2str(nfig)]);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for kk=1:nplo
                fprintf(fidTeX,'\\psfrag{%s}[1][][1][0]{%s}\n',deblank(vargin{indplot(l1)}),['$', deblank(texname{indplot(l1)}), '$']);
                l1=l1+1;
            end
            %         numberofvars=numberofvars-9;
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[scale=0.7] {%s} \n',ab);
            fprintf(fidTeX,'\\caption{Annualized Deterministic simulation %s:%s}', int2str(simul_nbr), int2str(nfig));
            fprintf(fidTeX,'\\label{Fig:DetSim:%s:%s}\n',simul_name,int2str(nfig));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
            fclose(fidTeX);
        end
        nplo=0;
    end
end

