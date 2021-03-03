function gr2perc_posterior(cname, lname, varargin)
global M_ oo_ options_

fname_ =M_.fname;
DirectoryName = CheckPath('Output');

iexo=[];
for i=1:length(varargin)
  iexo=[iexo, strmatch(varargin{i},M_.exo_names,'exact')];
end
iendo=[];
for j=1:length(cname),
  iendo=[iendo, strmatch(cname{j},M_.endo_names,'exact')];
end

nsubplo=9;
nn=3;

y=getIRFRuns(iendo,iexo);
for i=1:length(varargin)
  nfig=0;
  nplo=0;
  for j=1:length(cname),
    if length(iexo)==1 & length(iendo)==1;    
      y0=squeeze(y(:,:));
  elseif length(iexo)==1,
      y0=squeeze(y(:,j,:));
    elseif length(iendo)==1
      y0=squeeze(y(:,i,:));
    else
      y0=squeeze(y(:,j,i,:));
    end
    y1=cumsum(y0);
    for k=1:options_.irf,
        [MeanIRF(k,j,i),MedianIRF(k,j,i),VarIRF(k,j,i),HPDIRF(k,:,j,i),DistribIRF(k,:,j,i)] = ...
          posterior_moments(y1(k,:),0);
    end
    
    name = [lname{j} '_' varargin{i}];
    eval(['oo_.PosteriorIRF.dsge.Mean.' name ' = MeanIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.Median.' name ' = MedianIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.Var.' name ' = VarIRF(:,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.Distribution.' name ' = DistribIRF(:,:,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.HPDinf.' name ' = HPDIRF(:,1,j,i);']);
    eval(['oo_.PosteriorIRF.dsge.HPDsup.' name ' = HPDIRF(:,2,j,i);']);
    if max(y1) - min(y1) > 1e-10 ,
      nplo=nplo+1;
      if mod(nplo,nsubplo)==1,
        figure('name',['Orthogonalised shocks to ',varargin{i}, ' in perc''s'])
        nfig=nfig+1;
      end
      subplot(nn,nn,nplo)
      plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
      hold on
%       for k = 1:9
%         plot(1:options_.irf,DistribIRF(:,k,j,i),'-g','linewidth',0.5)
%       end
      for k = 1:2
        plot(1:options_.irf,HPDIRF(:,k,j,i),'-g','linewidth',0.5)
      end
      plot(1:options_.irf,MeanIRF(:,j,i),'-k','linewidth',1)
      xlim([1 options_.irf]);
      hold off
      title(lname{j},'interpreter','none')
    end
    if (mod(nplo,nsubplo)==0 | j==length(cname)) & nplo,
      saveas(gcf,[DirectoryName '\' fname_,'_Bayesian_IRF_',varargin{i},'_P',int2str(nfig)])
      eval(['print -depsc2 ' DirectoryName '\' fname_,'_Bayesian_IRF_',varargin{i},'_P',int2str(nfig)]);
      eval(['print -dpdf ' DirectoryName,'\',fname_,'_Bayesian_IRF_',varargin{i},'_P',int2str(nfig)]);
      close(gcf)
      nplo=0;
    end
  end
end
