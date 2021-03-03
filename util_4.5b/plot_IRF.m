function plot_IRF(xname, vargin, s1, texname, nr, nc, dirname, labels, Tend, init)
global M_ options_

fname_ =M_.fname;

if nargin<4 || isempty(texname) , texname=''; end
if nargin<7 || isempty(dirname), dirname='Output'; end
if nargin<8 || isempty(labels), xscale=1; end
if ~isfield(M_,'dname')
    M_.dname = M_.fname;
end
DirectoryName = CheckPath(dirname,M_.dname);
% DirectoryName =pwd;

if nargin<5 || isempty(nr)
    nr = 3;
end
if nargin<6 || isempty(nc)
    nc = 4;
end

nsub = nc*nr;
if nargin<9 || isempty(Tend)
    Tend = options_.irf;
else
    Tend = min(options_.irf, Tend);
end    

if nargin<10 || isempty(init), %#ok<*NOCOL>
    init = 1;
end

if init,
    delete([DirectoryName filesep,fname_ '_IRF.TeX']);
end

if nargin==0,
    return
end

if ~iscell(vargin)
    for j=1:size(vargin,1)
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

if ~iscell(xname)
    for j=1:size(xname,1);
        xname0{j}=deblank(xname(j,:));
    end
    xname=xname0;
end

iendo=[];
for i=1:length(vargin),
    itemp = strmatch(vargin{i},M_.endo_names,'exact');
    if isempty(itemp)
        error(['Variable ',vargin{i},' is not defined!'])
    else
        iendo=[iendo, itemp];
    end
end
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

if (nargin<4 && options_.TeX) || (isempty(texname) && options_.TeX)
    for j=1:length(vargin),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end

numberofvars=length(vargin);
for i=1:length(iexo),
nfig=0;
nplo=0;
indplot=[];
for j=1:length(vargin),
    
    name = [vargin{j} '_' deblank(M_.exo_names(iexo(i),:))];
    try
        MeanIRF=s1.(name);
    catch
        MeanIRF = 0;
    end
    if max(abs(MeanIRF)) > 1e-6 ,
        nplo=nplo+1;
        indplot=[indplot;j];
        if mod(nplo,nsub)==1,
            dyn_figure(options_.nodisplay,'name',['Orthogonalised shocks to ',deblank(M_.exo_names(iexo(i),:))]);
            nfig=nfig+1;
        end
        subplot(nr,nc,nplo)
        if nargin==8, 
            xscale=labels{j,3}; 
        end
        plot([1 Tend]/xscale,[0 0],'-r','linewidth',0.5);
        hold on,
        plot([1:Tend]/xscale,MeanIRF(1:Tend)*vscale(j)*exoscale(i),'-k','linewidth',2)
        xlim([1 Tend]/xscale);
        hold off
        if nargin==8, 
            text(0.82,-0.2,labels{j,1},'units','normalized');
            text(0.03,0.89,labels{j,2},'units','normalized');
        end
        if options_.TeX==0,
            title(texname{j},'interpreter','tex')
        else
            title(vargin{j},'interpreter','none')
        end
    end
    if (mod(nplo,nsub)==0 || j==length(vargin)) && nplo,
            dyn_saveas(gcf,[DirectoryName filesep fname_,'_IRF_',int2str(Tend),'_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
        if options_.TeX,
        texfilename=[DirectoryName filesep,fname_ '_IRF.TeX'];
            l1=(nc*nr)*(nfig-1)+1;
            
            %         texname0=texname(posi,:);
            %         fidTeX = fopen([M_.fname '_Fit_Plots.TeX'],'w+');
            %         files=dir('*_Fit*.eps');
            fidTeX = fopen([DirectoryName filesep,fname_ '_IRF.TeX'],'a+');
%         tempfilename=[DirectoryName filesep fname_,'_IRF_',deblank(M_.exo_names(iexo(i),:)),'_1'];
            %         tempfilename=deblank(tempfilename);
            %         [a,b,c]=fileparts(tempfilename);
            %         ab=strcat(a,'/',b);
        ab=([DirectoryName '/',fname_,'_IRF_',int2str(Tend),'_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for kk=1:nplo,
                fprintf(fidTeX,'\\psfrag{%s}[1][][1][0]{%s}\n',deblank(vargin{indplot(l1)}),['$', deblank(texname{indplot(l1)}), '$']);
                l1=l1+1;
            end
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
