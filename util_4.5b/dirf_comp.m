function dirf_comp(var_list_, s1, s2, vargin, texname)
global M_ options_
DirectoryName = CheckPath('Output');

fname_ = M_.fname;

if nargin<5 & options_.TeX,
for j=1:M_.endo_nbr,
  texname{j}=deblank(M_.endo_names_tex(j,:));
end
end

if ~iscell(vargin)
  for j=1:size(vargin,1);
    vargin0{j}=deblank(vargin(j,:));
  end
  vargin=vargin0;
end

iendo=[];
for i=1:size(var_list_,1)
  iendo=[iendo, strmatch(deblank(var_list_(i,:)),M_.endo_names,'exact')];
  if nargin<5 & options_.TeX,
    if strmatch('\log',texname{iendo(i)});
      texname{iendo(i)}=texname{iendo(i)}(2:end);
    end
  elseif options_.TeX,
    if strmatch('\log',texname{i});
      texname{i}=texname{i}(2:end);
    end
end
end


for i=1:length(vargin)
  nfigsav=0;
  nfig=0;
  nplo=0;
  for j=1:size(var_list_,1),
    if isfield(s1,[deblank(var_list_(j,:)),'_',vargin{i}]),
      nplo=nplo+1;
      if mod(nplo,9)==1,
        figure('name',['Comparison of orthogonalised shocks to ',vargin{i}])
        nfig=nfig+1;
      end
      subplot(3,3,nplo-9*(nfig-1))
      y=getfield(s1,[deblank(var_list_(j,:)),'_',vargin{i}]);
      if isfield(s2,[deblank(var_list_(j,:)),'_',vargin{i}])
        y=[y getfield(s2,[deblank(var_list_(j,:)),'_',vargin{i}])];
      else
        y=[y zeros(size(y))];
      end
      if ~exist('ys');
        ys=zeros([size(y), size(var_list_,1)]);
      end
      ys(:,:,j)=ys(:,:,j)+y;
      x=[1:length(y)];
      plot(x,y(:,1),'k',x,y(:,2),':k'), 
      if options_.TeX,
        if nargin<5
          title(deblank(texname{iendo(j)}),'interpreter','tex')
        else
          title(deblank(texname{j}),'interpreter','tex')
        end
      else    
        title(deblank(var_list_(j,:)),'interpreter','none')
      end
      x0=get(gca,'xlim');
      hold on, plot(x0, [0 0],'r')
    end
      if (mod(nplo,9)==0 | j==size(var_list_,1)) & nfigsav<nfig & nfig>0,
        nfigsav=nfig;
        saveas(gcf,[DirectoryName '\' fname_,'_IRF_comp_',vargin{i},int2str(nfig)])
        eval(['print -depsc2 ' DirectoryName '\' fname_,'_IRF_comp_',vargin{i},int2str(nfig)]);
        eval(['print -dpdf ' DirectoryName '\' fname_,'_IRF_comp_',vargin{i},int2str(nfig)]);
        close(gcf)
      end
  end
end


  nfigsav=0;
  nfig=0;
  nplo=0;
for j=1:size(var_list_,1),
      nplo=nplo+1;
      if mod(nplo,9)==1,
        figure('name',['Comparison of cumulated orthogonalised shocks'])
        nfig=nfig+1;
      end
      subplot(3,3,nplo-9*(nfig-1))
      plot(x,ys(:,1,j),'k',x,ys(:,2,j),':k'), title(deblank(var_list_(j,:)),'interpreter','none')
      x0=get(gca,'xlim');
      hold on, plot(x0, [0 0],'r')
      if (mod(nplo,9)==0 | j==size(var_list_,1)) & nfigsav<nfig & nfig>0,
        nfigsav=nfig;
        saveas(gcf,[DirectoryName '\' fname_,'_IRF_cum_',int2str(nfig)])
        eval(['print -depsc2 ' DirectoryName '\' fname_,'_IRF_cum_',int2str(nfig)]);
        eval(['print -dpdf ' DirectoryName '\' fname_,'_IRF_cum_',int2str(nfig)]);
        close(gcf)
      end
end

