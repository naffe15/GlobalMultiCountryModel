function plot_simul(vargin, texname, nr, nc, dirname, tlim, out, legtxt,simul_name,init)
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

iendo=[];
for i=1:length(vargin)
    itemp = strmatch(vargin{i},M_.endo_names,'exact');
    if isempty(itemp),
        error(['variable ',vargin{i},' does not exist']);
    else        
        iendo=[iendo, strmatch(vargin{i},M_.endo_names,'exact')];
    end
end
if (nargin<2) || (isempty(texname) )
    for j=1:length(vargin),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end
if init,
    if options_.TeX,
        delete([DirectoryName '\',fname_ '_simul.TeX']);
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

if nargin<6 || isempty(tlim),
    name = vargin{1};
    [z,zss]=dyn2vec(name);
    tlim=[1 length(z)];
end    

for j=1:length(vargin),
    
    name = vargin{j};
    [z,zss]=dyn2vec(name);
    if nargin >=7,
        yy=zeros(length(z),length(out));
        for jy=1:length(out),
            yy(:,jy)=getfield(out(jy),name);
        end        
        yy=yy(M_.maximum_exo_lag+1:end,:);
    end
    z=z(M_.maximum_exo_lag+1:end);
    if flag_scale,
        z=z-zss;
        zss=0;
    end
    nplo=nplo+1;
    indplot=[indplot;j];
    if mod(nplo,nsub)==1,
        dyn_figure(options_.nodisplay,'name',['Deterministic simulation: ', simul_name]);
        nfig=nfig+1;
    end
    subplot(nr,nc,nplo)
    plot(z*vscale(j),'-k','linewidth',1)
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
    plot(tlim,[zss zss],'r','linewidth',0.5);
    xlim(tlim); 
%     xtick = get(gca,'xtick');
%     xtick=xtick(find(mod(xtick,4)==0));
    xtick=1:4:tlim(2);
    xtick=xtick(xtick>=tlim(1));
    dtick = ceil(length(xtick)/4);
    ntick0 = ceil(dtick/2);
    xtick=xtick(ntick0:dtick:end);
%     if length(xtick)<3 
%         if tlim(end)>20
%         dxtick=floor(sqrt(tlim/4));
%         xtick=4:4:tlim(2);
%         xtick=xtick(5:5:end);
%         end
%     end
    set(gca,'xtick', xtick);
    set(gca,'xticklabel',(xtick-1)/4)
    xlabel('periods (yrs)')
    hold off
    if options_.TeX==0
        title(texname{j},'interpreter','tex')
    else
        title(vargin{j},'interpreter','none')
    end
    if (mod(nplo,nsub)==0 || j==length(vargin)) && nplo,
        dyn_saveas(gcf,[DirectoryName filesep fname_,'_simul_',simul_name,'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
        if options_.TeX
            texfilename=[DirectoryName filesep fname_ '_simul.TeX'];
            l1=(nc+nr)*(nfig-1)+1;
            
            fidTeX = fopen([DirectoryName filesep fname_ '_simul.TeX'],'a+');
            ab=([DirectoryName '/' fname_,'_simul_',simul_name,'_',int2str(nfig)]);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for kk=1:nplo
                fprintf(fidTeX,'\\psfrag{%s}[1][][1][0]{%s}\n',deblank(vargin{indplot(l1)}),['$', deblank(texname{indplot(l1)}), '$']);
                l1=l1+1;
            end
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
