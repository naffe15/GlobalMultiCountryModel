function myDefault

set(0,'DefaultFigurePaperType','a4')
%set(0,'DefaultFigurePaperOrientation','landscape')
set(0,'DefaultFigureColor',[0 0 0])
set(0,'DefaultTextColor',[1 1 1])
set(0,'DefaultAxesUnits','normalized')
set(0,'DefaultTextUnits','normalized')
set(0,'DefaultAxesLineStyleOrder', ...
  {'-',':','--','-.'})
set(0,'DefaultAxesColor',[0 0 0],...
    'DefaultAxesXColor', [1 1 1], ...
    'DefaultAxesYColor', [1 1 1], ...
    'DefaultAxesZColor', [1 1 1])
set(0, ...
      'DefaultAxesColorOrder', ...
      [1 1 0; ... %yellow
      1 1 1; ... % white
      0 1 1; ... % cyan
      1 0 1; ... % magenta
      1 0 0; ... % red
      0 1 0; ... % green
      0 0 1]) % blue
 

%set(0,'DefaultPatchmarkeredgecolor',[1 1 1])
set(0,'DefaultPatchEdgeColor',[1 1 1])
%map = colormap;
%colormap(map([end:-1:1],:));
%colormap(hsv);
