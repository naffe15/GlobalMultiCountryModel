function plot_conditionalIRF_post(xname, vargin, texname)
global M_ oo_ options_ bayestopt_

fname_ =M_.fname;
DirectoryName = CheckPath('Output');


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

y=getIRFRuns(iendo,iexo);
DirName = CheckPath('Metropolis');
load([DirName,'\',M_.fname,'_param_irf1.mat'])

for i=1:length(iexo),
  nfig=0;
  nplo=0;
  ipar=strmatch(xname{i},bayestopt_.name,'exact');
  for j=1:length(vargin),
    
    y0=squeeze(y(:,j,i,:));
    if ~isempty(ipar),
      for k=1:size(y0,1),
        x0=stand(stock(:,ipar));
        b=regress(y0(k,:)',[ones(500,1) x0]);
        y0(k,:)=y0(k,:)-x0'.*b(2);
      end
    end
    for k=1:options_.irf,
        [MeanIRF(k),MedianIRF(k),VarIRF(k),HPDIRF(k,:),DistribIRF(k,:)] = ...
          posterior_moments(y0(k,:),0);
    end
    name = [vargin{j} '_' deblank(M_.exo_names(iexo(i),:))];

    if max(abs(MeanIRF)) > 1e-6 ,
      nplo=nplo+1;
      if mod(nplo,9)==1,
        figure('name',['Orthogonalised shocks to ',deblank(M_.exo_names(iexo(i),:))])
        nfig=nfig+1;
      end
      HPDinf = HPDIRF(:,1);
      HPDsup = HPDIRF(:,2);
      subplot(3,3,nplo)
      patch([[1:options_.irf] [options_.irf:-1:1]],[HPDsup' HPDinf(end:-1:1)'],[0.75 0.75 0.75])
      hold on
      plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
%       for k = 1:9
%         plot(1:options_.irf,DistribIRF(:,k),'-g','linewidth',0.5)
%       end
      plot(1:options_.irf,MeanIRF,'-k','linewidth',1)
      xlim([1 options_.irf]);
      hold off
      if options_.TeX,
        title(texname{j},'interpreter','tex')
      else
        title(vargin{j},'interpreter','none')
      end
    end
    if (mod(nplo,9)==0 | j==length(vargin)) & nplo,
      saveas(gcf,[DirectoryName '\' fname_,'_Bayesian_cond_IRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)])
      eval(['print -depsc2 ' DirectoryName '\' fname_,'_Bayesian_cond_IRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
      eval(['print -dpdf ' DirectoryName,'\',fname_,'_Bayesian_cond_IRF_',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
      close(gcf)
      nplo=0;
    end
  end
end
