function plot_smooth_post(vname, nr, nc, texname, dirname)
global M_ oo_ options_

fname_ =M_.fname;
s1 = oo_.SmoothedVariables;

if nargin<5, dirname='Output'; end
DirectoryName = CheckPath(dirname,M_.dname);
if nargin<2,
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

if ~iscell(vname)
    for j=1:size(vname,1);
        vname0{j}=deblank(vname(j,:));
    end
    vname=vname0;
end

iendo=[];
for i=1:length(vname)
    iendo=[iendo, strmatch(vname{i},M_.endo_names,'exact')];
end

if (nargin<4) || (isempty(texname))
    for j=1:length(iendo),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end
delete([DirectoryName '\',fname_ '_PosteriorSmoothedVariables.TeX']);

% numberofvars=length(vargin);
nplo=0;
nfig=0;
for i=1:length(iendo),
    indplot=[];
%     for j=1:length(vargin),
%           Mean: [1x1 struct]
%           Median: [1x1 struct]
%              Var: [1x1 struct]
%     Distribution: [1x1 struct]
%           HPDinf: [1x1 struct]
%           HPDsup: [1x1 struct]

        name = [deblank(M_.endo_names(iendo(i),:))];
        eval(['Mean=s1.Mean.' name,';']);
        eval(['Distrib=s1.deciles.' name ';']);
        eval(['HPDsup=s1.HPDsup.' name ';']);
        eval(['HPDinf=s1.HPDinf.' name ';']);
        
        if max(abs(Mean)) > 1e-6 ,
            nplo=nplo+1;
            indplot=[indplot;1];
            if mod(nplo,nsub)==1,
                nfig=nfig+1;
                dyn_figure(options_.nodisplay,'name',['Posterior Smoothed Variables ',int2str(nfig)]);
            end
            subplot(nr,nc,nplo)
            patch([[1:options_.nobs] [options_.nobs:-1:1]],[HPDsup; HPDinf(end:-1:1)]*vscale(1),[0.75 0.75 0.75])
%             patch([[1:options_.nobs]' [1:options_.nobs]'],[HPDsup' HPDinf']*vscale(1),[0.75 0.75 0.75])
%             patch([HPDsup HPDinf]*vscale(1),[0.75 0.75 0.75])
            hold on
            plot([1 options_.nobs],[get_mean(name) get_mean(name)],'-r','linewidth',0.5);
%             plot([1 options_.nobs],[0 0],'-r','linewidth',0.5);
            %       for k = 1:9
            %         plot(1:options_.irf,DistribIRF(:,k),'-g','linewidth',0.5)
            %       end
            plot(1:options_.nobs,Mean(:)*vscale(1),'-k','linewidth',2)
            xlim([1 options_.nobs]);
            hold off
            if options_.TeX==0,
                title(texname{i},'interpreter','tex')
            else
                title(vname{i},'interpreter','none')
            end
        end
        if i==length(iendo) || ((mod(nplo,nsub)==0 ) && nplo),
            dyn_saveas(gcf,[DirectoryName '\' fname_,'_PosteriorSmoothedVariables_',int2str(nfig)],options_.nodisplay,options_.graph_format)
            if options_.TeX,
                texfilename=[DirectoryName '\',fname_ '_PosteriorSmoothedVariables.TeX'];
                l1=(nc+nr)*(nfig-1)+1;
                
                %         texname0=texname(posi,:);
                %         fidTeX = fopen([M_.fname '_Fit_Plots.TeX'],'w+');
                %         files=dir('*_Fit*.eps');
                fidTeX = fopen(texfilename,'a+');
                %         tempfilename=[DirectoryName '\' fname_,'_IRF_',deblank(M_.exo_names(iexo(i),:)),'_1'];
                %         tempfilename=deblank(tempfilename);
                %         [a,b,c]=fileparts(tempfilename);
                %         ab=strcat(a,'/',b);
                ab=([DirectoryName '/',fname_,'_PosteriorSmoothedVariables_',int2str(nfig)]);
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for kk=1:nplo,
                    fprintf(fidTeX,'\\psfrag{%s}[1][][1][0]{%s}\n',deblank(vname{indplot(l1)}),['$', deblank(texname{indplot(l1)}), '$']);
                    l1=l1+1;
                end
                %         numberofvars=numberofvars-9;
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,['\\includegraphics[scale=0.7] {%s} \n'],ab);
                fprintf(fidTeX,'\\caption{Posterior smoothed variables $%s$}',int2str(nfig));
                fprintf(fidTeX,'\\label{Fig:PosteriorSmoothedVariables:%s}\n',int2str(nfig));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,' \n');
                fclose(fidTeX);
            end
            nplo=0;
        end
    end
end
