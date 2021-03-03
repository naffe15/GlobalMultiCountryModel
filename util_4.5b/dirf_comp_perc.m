function dirf_comp_perc(cname, lname, s1, s2, vargin, texname)
global M_ options_

fname_=M_.fname;
DirectoryName = CheckPath('Output');

if nargin<6,
  for j=1:M_.endo_nbr,
    texname{j}=deblank(M_.endo_names_tex(j,:));
  end
end


for i=1:length(vargin)
  nfig=0;
  nplo=0;
  for j=1:length(cname),
    if isfield(s1,[cname{j},'_',vargin{i}]),
      nplo=nplo+1;
      if mod(nplo,9)==1,
        figure('name',['Perc. change comparison of orthogonalised shocks to ',vargin{i}])
        nfig=nfig+1;
      end
      subplot(3,3,nplo-9*(nfig-1))
      ytemp_=getfield(s1,[cname{j},'_',vargin{i}]);
      y=cumsum(ytemp_);
      ytemp_=getfield(s2,[cname{j},'_',vargin{i}]);
      y=[y cumsum(ytemp_)];
      if ~exist('ys');
        ys=zeros([size(y), length(cname)]);
      end
      ys(:,:,j)=ys(:,:,j)+y;
      x=[1:length(y)];
      plot(x,y(:,1),'k',x,y(:,2),':k'), 
      if options_.TeX,
        if nargin<6
          title(lname{j},'interpreter','none')
        else
          title(texname{j},'interpreter','tex')
        end
      else
        title(lname{j},'interpreter','none')
      end
      x0=get(gca,'xlim');
      hold on, plot(x0, [0 0],'r')
    end
    if (mod(nplo,9)==0 | j==length(cname)) & nfig>0,
      saveas(gcf,[DirectoryName '\' fname_,'_IRF_comp_perc_',vargin{i},int2str(nfig)])
      eval(['print -depsc2 ' DirectoryName '\' fname_,'_IRF_comp_perc_',vargin{i},int2str(nfig)]);
      eval(['print -dpdf ' DirectoryName '\' fname_,'_IRF_comp_perc_',vargin{i},int2str(nfig)]);
      close(gcf)
    end
  end
end

  nfigsav=0;
  nfig=0;
  nplo=0;
for j=1:length(cname),
      nplo=nplo+1;
      if mod(nplo,9)==1,
        figure('name',['Perc. comparison of cumulated orthogonalised shocks'])
        nfig=nfig+1;
      end
      subplot(3,3,nplo-9*(nfig-1))
      plot(x,ys(:,1,j),'k',x,ys(:,2,j),':k'), 
      if options_.TeX,
        title(texname{j},'interpreter','tex')
      else
        title(lname{j},'interpreter','none')
      end
      x0=get(gca,'xlim');
      hold on, plot(x0, [0 0],'r')
      if (mod(nplo,9)==0 | j==length(cname)) & nfigsav<nfig & nfig>0,
        nfigsav=nfig;
        saveas(gcf,[DirectoryName '\' fname_,'_IRF_cum_perc_',int2str(nfig)])
        eval(['print -depsc2 ' DirectoryName '\' fname_,'_IRF_cum_perc_',int2str(nfig)]);
        eval(['print -dpdf ' DirectoryName '\' fname_,'_IRF_cum_perc_',int2str(nfig)]);
        close(gcf)
      end
end

