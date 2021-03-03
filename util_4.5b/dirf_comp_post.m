function dirf_comp_post(var_list_, s1, s2, vargin, texname)
global M_ options_

fname_ = M_.fname;
DirectoryName = CheckPath('Output');

for i=1:length(vargin)
  nfigsav=0;
  nfig=0;
  nplo=0;
  for j=1:length(var_list_),
    y=getfield(s1.PosteriorIRF.dsge.Mean,[var_list_{j},'_',vargin{i}]);
    y=[y getfield(s2.PosteriorIRF.dsge.Mean,[var_list_{j},'_',vargin{i}])];
    if max(max(abs(y)))>1.e-10,
      nplo=nplo+1;
      if mod(nplo,9)==1,
        figure('name',['Comparison of orthogonalised shocks to ',vargin{i}])
        nfig=nfig+1;
      end
      subplot(3,3,nplo)
      yl=getfield(s1.PosteriorIRF.dsge.HPDinf,[var_list_{j},'_',vargin{i}]);
      yu=getfield(s1.PosteriorIRF.dsge.HPDsup,[var_list_{j},'_',vargin{i}]);
      
      yl=[yl getfield(s2.PosteriorIRF.dsge.HPDinf,[var_list_{j},'_',vargin{i}])];
      yu=[yu getfield(s2.PosteriorIRF.dsge.HPDsup,[var_list_{j},'_',vargin{i}])];
      x=[1:size(y,1)];
      patch([x, x(end:-1:1)], [yu(: ,1)', yl(end:-1:1,1)'],[0.75 0.75 0.75]); %,'facealpha',0.5);
      hold on,
      %patch([x, x(end:-1:1)], [yu(: ,2)', yl(end:-1:1,2)'],[0.75 0 0],'facealpha',0.5);
      plot(x,y(:,1),'k',x,y(:,2),'--k','LineWidth',1), 
      if options_.TeX & nargin==5
        title(texname{j},'interpreter','tex')
      else
        title(var_list_{j},'interpreter','none')
      end
      plot(x([1, end]), [0 0],'r')
      hold off,
    end
      if (mod(nplo,9)==0 | j==length(var_list_)) & nfigsav<nfig & nfig>0,
        nfigsav=nfig;
        saveas(gcf,[DirectoryName,'\',fname_,'_Bayesian_IRF_comp_',vargin{i},'_',int2str(nfig)])
        eval(['print -depsc2 ' DirectoryName,'\',fname_,'_Bayesian_IRF_comp_',vargin{i},'_',int2str(nfig)]);
        eval(['print -dpdf ' DirectoryName,'\',fname_,'_Bayesian_IRF_comp_',vargin{i},'_',int2str(nfig)]);
        close(gcf)
        nplo=0;
      end
  end
end

