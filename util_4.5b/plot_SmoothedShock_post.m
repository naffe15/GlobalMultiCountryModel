function plot_SmoothedShock_post(xname,  s1, texname, nr, nc, dirname)
global M_ oo_ options_

fname_ =M_.fname;

if nargin<7, dirname='Output'; end
DirectoryName = CheckPath(dirname,M_.dname);
if nargin<5,
    nr = 3;
    nc = 4;
end

nsub = nc*nr;


% if ~iscell(vargin)
%     for j=1:size(vargin,1);
%         vargin0{j,1}=deblank(vargin(j,:));
%     end
%     vargin=vargin0;
% end
% if size(vargin,2)==2,
%     temp1=vargin(:,2);
%     vargin=vargin(:,1);
%     for j=1:length(vargin),
%         vscale(j,1)=temp1{j,1};
%     end
% else
    vscale=1;
%     
% end

if ~iscell(xname)
    for j=1:size(xname,1);
        xname0{j}=deblank(xname(j,:));
    end
    xname=xname0;
end

iexo=[];
for i=1:length(xname)
    iexo=[iexo, strmatch(xname{i},M_.exo_names,'exact')];
end

if (nargin<3) || (isempty(texname))
    for j=1:length(iexo),
        texname{j}=deblank(M_.exo_names_tex(iexo(j),:));
    end
end
delete([DirectoryName '\',fname_ '_PosteriorSmoothedShocks.TeX']);

% numberofvars=length(vargin);
for i=1:length(iexo),
    nfig=0;
    nplo=0;
    indplot=[];
%     for j=1:length(vargin),
%           Mean: [1x1 struct]
%           Median: [1x1 struct]
%              Var: [1x1 struct]
%     Distribution: [1x1 struct]
%           HPDinf: [1x1 struct]
%           HPDsup: [1x1 struct]

        name = [deblank(M_.exo_names(iexo(i),:))];
        eval(['MeanSmoothedShocks=s1.Mean.' name,';']);
        eval(['DistribSmoothedShocks=s1.deciles.' name ';']);
        eval(['HPDsup=s1.HPDsup.' name ';']);
        eval(['HPDinf=s1.HPDinf.' name ';']);
        
        if max(abs(MeanSmoothedShocks)) > 1e-6 ,
            nplo=nplo+1;
            indplot=[indplot;1];
            if mod(nplo,nsub)==1,
                nfig=nfig+1;
                dyn_figure(options_.nodisplay,'name',['Posterior Smoothed Shocks ',int2str(nfig)]);
            end
            subplot(nr,nc,nplo)
            patch([[1:options_.nobs]' [options_.nobs:-1:1]'],[HPDsup' HPDinf(end:-1:1)']*vscale(1),[0.75 0.75 0.75])
%             patch([[1:options_.nobs]' [1:options_.nobs]'],[HPDsup' HPDinf']*vscale(1),[0.75 0.75 0.75])
%             patch([HPDsup HPDinf]*vscale(1),[0.75 0.75 0.75])
            hold on
            plot([1 options_.nobs],[0 0],'-r','linewidth',0.5);
            %       for k = 1:9
            %         plot(1:options_.irf,DistribIRF(:,k),'-g','linewidth',0.5)
            %       end
            plot(1:options_.nobs,MeanSmoothedShocks(:)*vscale(1),'-k','linewidth',1)
            xlim([1 options_.nobs]);
            hold off
            if options_.TeX==0,
                title(texname{i},'interpreter','tex')
            else
                title(texname{i},'interpreter','none')
            end
        end
        if (mod(nplo,nsub)==0 ) && nplo,
            dyn_saveas(gcf,[DirectoryName '\' fname_,'_PosteriorSmoothedShocks_',int2str(nfig)],options_.nodisplay,options_.graph_format)
            if options_.TeX,
                texfilename=[DirectoryName '\',fname_ '_PosteriorSmoothedShocks.TeX'];
                l1=(nc+nr)*(nfig-1)+1;
                
                %         texname0=texname(posi,:);
                %         fidTeX = fopen([M_.fname '_Fit_Plots.TeX'],'w+');
                %         files=dir('*_Fit*.eps');
                fidTeX = fopen(texfilename,'a+');
                %         tempfilename=[DirectoryName '\' fname_,'_IRF_',deblank(M_.exo_names(iexo(i),:)),'_1'];
                %         tempfilename=deblank(tempfilename);
                %         [a,b,c]=fileparts(tempfilename);
                %         ab=strcat(a,'/',b);
                ab=([DirectoryName '/',fname_,'_PosteriorSmoothedShocks_',int2str(nfig)]);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for kk=1:nplo,
                    fprintf(fidTeX,'\\psfrag{%s}[1][][1][0]{%s}\n',deblank(xname{indplot(l1)}),['$', deblank(texname{indplot(l1)}), '$']);
                    l1=l1+1;
                end
                %         numberofvars=numberofvars-9;
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,['\\includegraphics[scale=0.7] {%s} \n'],ab);
                fprintf(fidTeX,'\\caption{Posterior smoothed shocks $%s$}',int2str(nfig));
                fprintf(fidTeX,'\\label{Fig: Posterior Smoothed Shocks Plot:%s}\n',int2str(nfig));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,' \n');
                fclose(fidTeX);
            end
            nplo=0;
        end
    end
end
