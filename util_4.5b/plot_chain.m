function plot_chain
global M_ options_ estim_params_ bayestopt_

npar = estim_params_.np+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nvx;

nblck = options_.mh_nblck; 
DirectoryName = CheckPath('metropolis',M_.dname);
OutDirectoryName = CheckPath('output',M_.dname);
a=dir([DirectoryName filesep  M_.fname '_mh_history*']);
load([DirectoryName filesep  M_.fname '_mh_history_' int2str(length(a)-1)]);
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));

ifig=0;
for j=0:npar
  if mod(j+1,9)==1
    ifig=ifig+1;
    iplot=0;
    hh=dyn_figure(options_.nodisplay,'Name',['Markov Chain plot ',int2str(ifig)]);
  end
  iplot=iplot+1;
  Draws = GetAllPosteriorDraws(j,1,1,TotalNumberOfMhFiles,TotalNumberOfMhDraws);
  Draws=reshape(Draws,TotalNumberOfMhDraws,nblck);
  if j>0
      dtmp = Draws(ceil(TotalNumberOfMhDraws*options_.mh_drop):end,:);
      m1 = mean(dtmp(:));
      sd = std(dtmp(:));
  %     Draws = (Draws-m1)./sd;
  end
  subplot(3,3,iplot)
  plot(Draws,'linewidth',2)
  if j>0
      co = get(gca,'ColorOrder');
      set(gca, 'ColorOrder', co(1:4,:));
      hold on, plot(cumsum(Draws)./kron([1:length(Draws)]',ones(1,size(Draws,2))),':') %,'--r')
      plot([0 size(Draws,1)],[m1 m1],'--k')
      title(bayestopt_.name{j},'interpreter','none')
      hold off
  else
      title('log-posterior kernel','interpreter','none')
  end
  if mod(j+1,9)==0 || j==npar
    dyn_saveas(hh,[OutDirectoryName filesep M_.fname '_chain_' int2str(ifig)],options_.nodisplay,options_.graph_format)
  end
  
end


