function plot_IRF_post(xname, vargin, s1, texname, nr, nc, dirname)
global M_ oo_ options_

fname_ =M_.fname;

if nargin<7, dirname='Output'; end
DirectoryName = CheckPath(dirname,M_.dname);
if nargin<5,
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
for i=1:length(vargin)
    iendo=[iendo, strmatch(vargin{i},M_.endo_names,'exact')];
end
if nargin<3,
    for j=1:length(iendo),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end

iexo=[];
for i=1:length(xname)
    iexo=[iexo, strmatch(xname{i},M_.exo_names,'exact')];
end

if (nargin<4 && options_.TeX) || (isempty(texname) && options_.TeX)
    for j=1:length(vargin),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end
if ~ismac
    delete([DirectoryName '\',fname_ '_PosteriorIRF.TeX']);
else
    delete([DirectoryName '/',fname_ '_PosteriorIRF.TeX']);
end
numberofvars=length(vargin);
for i=1:length(iexo),
    nfig=0;
    nplo=0;
    indplot=[];
    for j=1:length(vargin),
        
        name = [vargin{j} '_' deblank(M_.exo_names(iexo(i),:))];
        eval(['MeanIRF=s1.Mean.' name,';']);
        eval(['DistribIRF=s1.deciles.' name ';']);
        eval(['HPDsup=s1.HPDsup.' name ';']);
        eval(['HPDinf=s1.HPDinf.' name ';']);
        irf_length = length(MeanIRF);
        
        if max(abs(MeanIRF)) > 1e-6 ,
            nplo=nplo+1;
            indplot=[indplot;j];
            if mod(nplo,nsub)==1,
                dyn_figure(options_.nodisplay,'name',['Orthogonalised shocks to ',deblank(M_.exo_names(iexo(i),:))]);
                nfig=nfig+1;
            end
            subplot(nr,nc,nplo)
            patch([[1:irf_length] [irf_length:-1:1]],[HPDsup; HPDinf(end:-1:1)]'*vscale(j),[0.75 0.75 0.75])
            hold on
            plot([1 irf_length],[0 0],'-r','linewidth',0.5);
            %       for k = 1:9
            %         plot(1:irf_length,DistribIRF(:,k),'-g','linewidth',0.5)
            %       end
            plot(1:irf_length,MeanIRF(:)*vscale(j),'-k','linewidth',1)
            xlim([1 irf_length]);
            hold off
            if options_.TeX==0,
                title(texname{j},'interpreter','tex')
            else
                title(vargin{j},'interpreter','none')
            end
        end
        if (mod(nplo,nsub)==0 || j==length(vargin)) && nplo,
            if ~ismac
                dyn_saveas(gcf,[DirectoryName '\' fname_,'_PosteriorIRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
            else
                dyn_saveas(gcf,[DirectoryName '/' fname_,'_PosteriorIRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
            end
                if options_.TeX,
                texfilename=[DirectoryName '\',fname_ '_PosteriorIRF.TeX'];
                l1=(nc+nr)*(nfig-1)+1;
                
                %         texname0=texname(posi,:);
                %         fidTeX = fopen([M_.fname '_Fit_Plots.TeX'],'w+');
                %         files=dir('*_Fit*.eps');
                fidTeX = fopen(texfilename,'a+');
                %         tempfilename=[DirectoryName '\' fname_,'_IRF_',deblank(M_.exo_names(iexo(i),:)),'_1'];
                %         tempfilename=deblank(tempfilename);
                %         [a,b,c]=fileparts(tempfilename);
                %         ab=strcat(a,'/',b);
                ab=([DirectoryName '/',fname_,'_PosteriorIRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for kk=1:nplo,
                    fprintf(fidTeX,'\\psfrag{%s}[1][][1][0]{%s}\n',deblank(vargin{indplot(l1)}),['$', deblank(texname{indplot(l1)}), '$']);
                    l1=l1+1;
                end
                %         numberofvars=numberofvars-9;
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,['\\includegraphics[scale=0.7] {%s} \n'],ab);
                fprintf(fidTeX,'\\caption{Impulse responses to a shock to $%s$}',[deblank(M_.exo_names_tex(iexo(i),:))]);
                fprintf(fidTeX,'\\label{Fig: Posterior IRF Plot:%s}\n',int2str(j));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,' \n');
                fclose(fidTeX);
            end
            nplo=0;
        end
    end
end
