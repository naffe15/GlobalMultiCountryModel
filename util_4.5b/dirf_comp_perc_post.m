function dirf_comp_perc_post(lname, s1, s2, varargin)
global M_ options_

fname_=M_.fname;
DirectoryName = CheckPath('Output');
in=size(lname);
if min(in)==2,
  texname=lname(2,:);
  lname=lname(1,:);
else
  for j=1:length(lname),
    iendo=strmatch(lname{j},M_.endo_names,'exact');
    texname{j}=deblank(M_.endo_names_tex(iendo,:));
  end
end

for i=1:length(varargin)
  nfig=0;
  nplo=0;
  for j=1:length(lname),
    y=getfield(s1.PosteriorIRF.Mean,[lname{j},'_',varargin{i}]);
    y=[y getfield(s2.PosteriorIRF.Mean,[lname{j},'_',varargin{i}])];
    if max(max(abs(y)))>1e-10 & max(max(y)-min(y))>1e-10;
      nplo=nplo+1;
      if mod(nplo,9)==1,
        figure('name',['Perc. change comparison of orthogonalised shocks to ',varargin{i}])
        nfig=nfig+1;
      end
      subplot(3,3,nplo)
      yl=getfield(s1.PosteriorIRF.HPDinf,[lname{j},'_',varargin{i}]);
      yl=[yl getfield(s2.PosteriorIRF.HPDinf,[lname{j},'_',varargin{i}])];
      yu=getfield(s1.PosteriorIRF.HPDsup,[lname{j},'_',varargin{i}]);
      yu=[yu getfield(s2.PosteriorIRF.HPDsup,[lname{j},'_',varargin{i}])];
      x=[1:length(y)];
%       patch([x, x(end:-1:1)], [yu(: ,1)', yl(end:-1:1,1)'],[0 0.75 0],'facealpha',0.5);
%       hold on,
%       patch([x, x(end:-1:1)], [yu(: ,2)', yl(end:-1:1,2)'],[0.75 0 0],'facealpha',0.5);
      plot(x,y(:,1),'k',x,y(:,2),':k'), 
      hold on,
      if options_.TeX,
        title(texname{j},'interpreter','tex')        
      else
        title(lname{j},'interpreter','none')
      end
      plot(x([1,end]), [0 0],'r')
      hold off,
    end
    if (mod(nplo,9)==0 | j==length(lname)) & nfig>0,
      saveas(gcf,[DirectoryName,'/',fname_,'_Bayesian_IRF_comp_perc_',varargin{i},'_',int2str(nfig)])
      eval(['print -depsc2 ' DirectoryName,'/',fname_,'_Bayesian_IRF_comp_perc_',varargin{i},'_',int2str(nfig)]);
      eval(['print -dpdf ' DirectoryName,'/',fname_,'_Bayesian_IRF_comp_perc_',varargin{i},'_',int2str(nfig)]);
      close(gcf)
      nplo=0;
    end
  end
end

