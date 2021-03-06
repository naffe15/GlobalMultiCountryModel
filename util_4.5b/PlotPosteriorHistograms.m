function PlotPosteriorHistograms()
% stephane.adjemian@ens.fr [09-09-2005]
global estim_params_ M_ options_ bayestopt_ oo_

OutputDirectoryName = CheckPath('Output');

TeX   	= options_.TeX;
nblck 	= options_.mh_nblck;
nvx   	= estim_params_.nvx;
nvn   	= estim_params_.nvn;
ncx   	= estim_params_.ncx;
ncn   	= estim_params_.ncn;
np    	= estim_params_.np ;
npar   	= nvx+nvn+ncx+ncn+np;

MaxNumberOfPlotPerFigure = 9;% The square root must be an integer!
nn = sqrt(MaxNumberOfPlotPerFigure);

figurename = 'Priors and posteriors hitograms';

if TeX    
  fidTeX = fopen([OutputDirectoryName '/' M_.fname '_PriorsAndPosteriorsHist.TeX'],'w');
  fprintf(fidTeX,'%% TeX eps-loader file generated by PlotPosteriorHistograms.m (Dynare).\n');
  fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
  fprintf(fidTeX,' \n');
end

figunumber = 0;
subplotnum = 0;

for i=1:npar
  subplotnum = subplotnum+1;
  if subplotnum == 1
    figunumber = figunumber+1;
    if options_.nodisplay
      hfig = figure('Name',figurename,'Visible','off');
    else
      hfig = figure('Name',figurename);
    end
  end
  if subplotnum == 1
    if TeX
      TeXNAMES = [];
    end
    NAMES = [];
  end
  [nam,texnam] = get_the_name(i,TeX);
  NAMES = strvcat(NAMES,nam);
  if TeX
    TeXNAMES = strvcat(TeXNAMES,texnam);
  end
  hist
  %[x2,f2,abscissa,dens,binf2,bsup2] = draw_prior_density(i);
  top2 = max(f2); 
  top1 = max(f1);
  top0 = max([top1;top2]);
  binf1 = x1(1);
  bsup1 = x1(end);
  borneinf = min(binf1,binf2);
  bornesup = max(bsup1,bsup2);
  subplot(nn,nn,subplotnum)
  hh = plot(x2,f2,'-k','linewidth',2);
  set(hh,'color',[0.7 0.7 0.7]);
  hold on;
  plot(x1,f1,'-k','linewidth',2);
  plot( [pmode pmode], [0.0 1.1*top0], '--g', 'linewidth', 2);
  box on;
  axis([borneinf bornesup 0 1.1*top0]);
  title(nam,'Interpreter','none');
  hold off;
  drawnow
  if subplotnum == MaxNumberOfPlotPerFigure | i == npar;
    eval(['print -depsc2 ' OutputDirectoryName '/' M_.fname '_PriorsAndPosteriors' int2str(figunumber)]);
    eval(['print -dpdf ' OutputDirectoryName '/' M_.fname '_PriorsAndPosteriors' int2str(figunumber)]);
    if options_.nodisplay, 
      set(hfig,'Visible','on');
    end
    saveas(hfig,[OutputDirectoryName '/' M_.fname '_PriorsAndPosteriors' int2str(figunumber) '.fig']);
    if TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for j = 1:size(NAMES,1)
	fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(j,:)),deblank(TeXNAMES(j,:)));
      end    
      fprintf(fidTeX,'\\centering\n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_PriorsAndPosteriors%s}\n',M_.fname,int2str(figunumber));
      fprintf(fidTeX,'\\caption{Priors and posteriors.}');
      fprintf(fidTeX,'\\label{Fig:PriorsAndPosteriors:%s}\n',int2str(figunumber));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,' \n');
      if i == npar
	fprintf(fidTeX,'%% End of TeX file.\n');
	fclose(fidTeX);
      end
    end
    if options_.nodisplay, 
      close(hfig), 
    end
    subplotnum = 0;
  end
end