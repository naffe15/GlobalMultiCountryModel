function PosteriorSmoother(type)
% stephane.adjemian@ens.fr [09-25-2005]
global options_ estim_params_ oo_ M_

nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
naK = length(options_.filter_step_ahead);
%%
endo_nbr=M_.endo_nbr;
exo_nbr=M_.exo_nbr;
nvobs     = size(options_.varobs,1);
%%
CheckPath('Plots/');
DirectoryName = CheckPath('metropolis');
%%
%%
varlist = options_.varlist;
if isempty(varlist)
  varlist = M_.endo_names;
  SelecVariables = transpose(1:M_.endo_nbr);
  nvar = M_.endo_nbr;
else
  nvar = size(varlist,1);
  SelecVariables = [];
  for i=1:nvar
    if ~isempty(strmatch(varlist(i,:),M_.endo_names,'exact'))
      SelecVariables = [SelecVariables;strmatch(varlist(i,:),M_.endo_names,'exact')];
    end
  end
end

h = waitbar(0,'Load Bayesian smoother...');
filsmooth = dir([DirectoryName '/' M_.fname '_smooth*.mat']);
y0=[];
for j=1:length(filsmooth),
  load([DirectoryName '/' M_.fname '_smooth',num2str(j),'.mat']);
  if isempty(y0),
    y0=stock;
  else
    %y0(:,:,size(y0,3):size(y0,3)+size(stock,3))=stock;
    y0=cat(3,y0,stock);
  end
end
if nvx,
  filinno = dir([DirectoryName '/' M_.fname '_inno*.mat']);
end
if nvn
  filerror = dir([DirectoryName '/' M_.fname '_error*.mat']);
end

if naK
  filfilt = dir([DirectoryName '/' M_.fname '_filter*.mat']);
end

filparam = dir([DirectoryName '/' M_.fname '_param*.mat']);

for b=1:B
  %deep = GetOneDraw(NumberOfDraws,FirstMhFile,LastMhFile,FirstLine,MAX_nruns,DirectoryName);
  deep = GetOneDraw(type);
  if irun1 > MAX_nsmoo | b == B
    if b == B
      stock_smooth = stock_smooth(:,:,1:irun1-1);
    end
    stock = stock_smooth;
    save([DirectoryName '/' M_.fname '_smooth' int2str(ifil1)],'stock');
    ifil1 = ifil1 + 1;
    irun1 = 1;
  end
  
  if nvx & (irun2 > MAX_ninno | b == B)
    if b == B
      stock_innov = stock_innov(:,:,1:irun2-1);
    end
    stock = stock_innov;
    save([DirectoryName '/' M_.fname '_inno' int2str(ifil2)],'stock');
    ifil2 = ifil2 + 1;
    irun2 = 1;
  end
  
  if nvn & (irun3 > MAX_error | b == B)
    if b == B
      stock_error = stock_error(:,:,1:irun3-1);
    end
    stock = stock_error;
    save([DirectoryName '/' M_.fname '_error' int2str(ifil3)],'stock');
    ifil3 = ifil3 + 1;
    irun3 = 1;
  end
  
  if naK & (irun4 > MAX_naK | b == B)
    if b == B
      stock_filter = stock_filter(:,:,:,1:irun4-1);
    end
    stock = stock_filter;
    save([DirectoryName '/' M_.fname '_filter' int2str(ifil4)],'stock');
    ifil4 = ifil4 + 1;
    irun4 = 1;
  end
  
  if irun5 > MAX_nruns | b == B
    if b == B
      stock_param = stock_param(1:irun5-1,:);
    end
    stock = stock_param;
    save([DirectoryName '/' M_.fname '_param' int2str(ifil5)],'stock');
    ifil5 = ifil5 + 1;
    irun5 = 1;
  end
  
  waitbar(b/B,h);
end
close(h)

