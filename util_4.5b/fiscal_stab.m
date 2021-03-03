function fiscal_stab(th_std_mh, vname, tname, aname)
global M_

clear fisc
if nargin<4, aname=''; end

fname_=M_.fname;

nv = size(th_std_mh,1);


for j=2:size(th_std_mh,3),
  figure
  fisc=squeeze((th_std_mh(:,:,j)-th_std_mh(:,:,1))./th_std_mh(:,:,1));
  fisc=fisc';
  fiscpol(:,j-1)=mean(fisc)';
  fiscpol_std(:,j-1)=std(fisc)';
  fisc=fisc.*100;
  fy = fisc;
  st=std(fy);
  ij=find(st~=0);
  %boxplot(fisc,'label',{' ',' ',' ',' ',' ',' '})
  boxplot(fy(:,ij),1,'.',[],3)
  hold on,
  plot([1 length(ij)], [0 0], 'g')
  yl=get(gca,'ylim');
  for i=1:length(ij),
    %text(i,(-1)^i*max((-1)^i*fisc(:,i)),cname{i})
    if abs(min(fisc(:,ij(i)))-yl(1))>abs(max(fisc(:,ij(i)))-yl(2)),
      text(i,min(fisc(:,ij(i))),[vname{ij(i)},' '],'rotation',90,'verticalalignment','middle','horizontalalignment','right')
    else
      text(i,max(fisc(:,ij(i))),[' ',vname{ij(i)}],'rotation',90,'verticalalignment','middle','horizontalalignment','left')
    end
  end
  title(tname{j-1})
  xlabel('')
  ylabel('% stabilisation effect')
  set(gca,'xticklabel','')
%   saveas(gcf,[fname_,'_fiscal',aname,num2str(j-1)])
%   eval(['print -depsc2 ' fname_,'_fiscal',aname,num2str(j-1)]);
%   eval(['print -dpdf ' fname_,'_fiscal',aname,num2str(j-1)]);
  disp(tname{j-1})
  sprintf('%5.2f%% (%5.2f%%) \n',[fiscpol(:,j-1)'; fiscpol_std(:,j-1)'].*100)
end

aa=squeeze(mean(th_std_mh,2));
save([fname_,'_fiscal_mh',aname], 'aa', 'fiscpol', 'fiscpol_std', '-append')

