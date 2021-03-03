function gr2perc(cname, lname, varargin)
global M_ options_

fname_ =M_.fname;
DirectoryName = CheckPath('Outputxx');

for i=1:length(varargin)
  nfig=0;
  nplo=0;
  for j=1:length(cname),
    if evalin('base',['exist(''',cname{j},'_',varargin{i},''')']),
      nplo=nplo+1;
      ytemp_ = evalin('base',[cname{j},'_',varargin{i}]);
      %y=exp(cumsum(ytemp_))-1;
      y=cumsum(ytemp_);
      if mod(nplo,9)==1,
        figure('name',['Orthogonalised shocks to ',varargin{i}, ' in perc''s'])
        nfig=nfig+1;
      end
      subplot(3,3,nplo-9*(nfig-1))
      plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
      hold on,
      plot([1:options_.irf], y,'k'), 
      assignin('base',[lname{j},'_',varargin{i}],y);
      title(lname{j},'interpreter','none')
      hold off,
    end
    if mod(nplo,9)==0 | j==length(cname),
      saveas(gcf,[DirectoryName '\' fname_,'_IRF_',varargin{i},'_P',int2str(nfig)])
      eval(['print -depsc2 ' DirectoryName '\' fname_,'_IRF_',varargin{i},'_P',int2str(nfig)]);
      eval(['print -dpdf ' DirectoryName '\' fname_,'_IRF_',varargin{i},'_P',int2str(nfig)]);
      close(gcf)
    end
  end
end
