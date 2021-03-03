function print_params_shocks_xls(fname,sheet,titcell);

global M_ bayestopt_

if isempty(bayestopt_)
  pnam=[];
else
  pnam=bayestopt_.name;
end

if nargin<3,
    titcell = pwd;
end

if ismac
   fix_xlsread_MACOS(fname); 
end
[NUMERIC,TXT,RAW]=xlsread(fname,sheet);

nraw=size(TXT,1);
icol = strmatch('ESTIMATED',TXT(1,:)');
icol1 = icol(1);
icol = icol(end);
lnam = [{'';'';''};TXT(4:end,icol)];
RAW = [RAW(1:nraw,1:icol-1) RAW(1:nraw,icol-2:icol-1) RAW(1:nraw,icol:end)];
for j=1:size(RAW,1),
  RAW{j,icol}='';
  RAW{j,icol+1}=[];
end


sd=sqrt(diag(M_.Sigma_e));

iraw = 3;
RAW{iraw,icol+1}  = titcell;
for j=1:M_.exo_nbr,
    
  jraw = iraw+1;
 
  iraw = strmatch(deblank(M_.exo_names(j,:)),lnam,'exact');
  if isempty(iraw),
     RAW = [RAW(1:jraw,:);RAW(jraw:end,:)] ;
     for jc = 1:size(RAW,2),
         RAW{jraw,jc}='';
     end
     RAW{jraw,icol+2}=deblank(M_.exo_names(j,:));
     RAW{jraw,icol1}=deblank(M_.exo_names(j,:));
     lnam = [{'';'';''};RAW(4:end,icol+2)];
     iraw = jraw;
  end
  
  
 RAW{iraw,icol+1}=sd(j);
  if strmatch(deblank(M_.exo_names(j,:)),pnam,'exact'),
    RAW{iraw,icol}='**';
    
  end
end

for j=1:length(M_.params)
  jraw = iraw + 1;
 
  iraw = strmatch(deblank(M_.param_names(j,:)),lnam,'exact');
  if isempty(iraw),
      if size(RAW,1)>=jraw,
     RAW = [RAW(1:jraw,:);RAW(jraw:end,:)] ;
      else
     RAW = [RAW;RAW(end,:)] ;
      end
     for jc = 1:size(RAW,2),
         RAW{jraw,jc}='';
     end
     RAW{jraw,icol+2}=deblank(M_.param_names(j,:));
     RAW{jraw,icol1}=deblank(M_.param_names(j,:));
     lnam = [{'';'';''};;RAW(4:end,icol+2)];
     iraw = jraw;
  end
  
  RAW{iraw,icol+1}=M_.params(j);
  if strmatch(deblank(M_.param_names(j,:)),pnam,'exact'),
    RAW{iraw,icol}='**';
  end
end
if ~ismac
    [SUCCESS,MESSAGE]=xlswrite(fname,RAW,sheet);
else
    [SUCCESS]=xlswrite_MACOS(fname,RAW,sheet);
end

