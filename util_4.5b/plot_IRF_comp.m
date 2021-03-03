function plot_IRF_comp(xname, vargin, s1, texname, nr, nc, dirname)
global M_ oo_ options_

if nargin<7 || isempty(dirname), dirname='Output'; end
fname_ =M_.fname;
if ~isfield(M_,'dname'),
    dname_=M_.fname;
else
    dname_=M_.dname;
end
DirectoryName = CheckPath(dirname,dname_);

if nargin<5,
    nr = 3;
    nc = 4;
end

nsub = nc*nr;

 
if ~iscell(vargin)
  for j=1:size(vargin,1);
    vargin0{j}=deblank(vargin(j,:));
  end
  vargin=vargin0;
end
if ~iscell(xname)
  for j=1:size(xname,1);
    xname0{j}=deblank(xname(j,:));
  end
  xname=xname0;
end

iendo=[];
for i=1:length(vargin)
  iendo=[iendo, strmatch(vargin{i},M_.endo_names,'exact')];
end
iexo=[];
for i=1:length(xname)
  iexo=[iexo, strmatch(xname{i},M_.exo_names,'exact')];
end

if (nargin<4 && options_.TeX) || (isempty(texname) && options_.TeX)
    for j=1:length(vargin),
        texname{j}=deblank(M_.endo_names_tex(iendo(j),:));
    end
end


for i=1:length(iexo),
  nfig=0;
  nplo=0;
  for j=1:length(vargin),
    clear MeanIRF,
    name = [vargin{j} '_' deblank(M_.exo_names(iexo(i),:))];
        for jcomp=1:length(s1)
    try
    eval(['MeanIRF(jcomp,:)=s1{jcomp}.' name,';']);
    catch
      MeanIRF(jcomp,:) = NaN;
    end
        end
    if max(max(abs(MeanIRF))) > 1e-6 ,
      irf_=size(MeanIRF,2);
      nplo=nplo+1;
      if mod(nplo,nsub)==1,
        dyn_figure(options_.nodisplay,'name',['Orthogonalised shocks to ',deblank(M_.exo_names(iexo(i),:))]);
        nfig=nfig+1;
      end
      subplot(nr,nc,nplo)
      plot([1 irf_],[0 0],'-r','linewidth',0.5);
      hold on,
      h=plot(1:irf_,MeanIRF','linewidth',1);
      set(h(1),'color','k')
      set(h(2),'color','k','linestyle',':')
      if length(h)==3,
          set(h(3),'color','k','linestyle','--')
      end
      xlim([1 irf_]);
      hold off
      if options_.TeX==0,
          title(texname{j},'interpreter','tex')
      else
          title(vargin{j},'interpreter','none')
      end
    end
    if (mod(nplo,nsub)==0 | j==length(vargin)) & nplo,
      dyn_saveas(gcf,[DirectoryName '\' fname_,'_IRF_comp',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)],options_.nodisplay,options_.graph_format)
%       eval(['print -depsc2 ' DirectoryName '\' fname_,'_IRF_comp',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
%       eval(['print -dpdf ' DirectoryName,'\',fname_,'_IRF_comp',deblank(M_.exo_names(iexo(i),:)),'_',int2str(nfig)]);
%       close(gcf)
      nplo=0;
    end
  end
end
