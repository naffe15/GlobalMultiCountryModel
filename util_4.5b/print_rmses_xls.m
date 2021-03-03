function print_rmses_xls(rmse,r2,lnam,rmseK,lnamK,fname,sheet,titcell);

global oo_

if nargin<8,
    titcell = pwd;
end

if ismac
   fix_xlsread_MACOS(fname); 
end

[NUMERIC,TXT,RAW]=xlsread(fname,sheet);
% if ismac
%    temp = cell(2, size(RAW,2));
%    RAW =[temp; RAW]; 
% end    
iraw = strmatch('log-dens',TXT);
iraw=iraw(1);
TXT=TXT(iraw-1:end,:);
icol = strmatch('log-dens',TXT(2,:)');
icol1 = icol(1);
icol = icol(end);
iraw = strmatch('1-step',TXT(:,2));
irawK = strmatch('4-step',TXT(:,2));
irawr2 = strmatch('r2',TXT(:,2));
vnam = TXT(iraw+2:irawK-1,icol);
vnamK = TXT(irawK+1:irawr2-1,icol);
RAW = [RAW(:,1:icol-1) RAW(:,icol-1:end)];
for j=1:size(RAW,1),
  RAW{j,icol}=[];
end


RAW{3,icol}  = titcell;
RAW{irawK+2,icol}  = titcell;
try
RAW{4,icol}  = oo_.MarginalDensity.LaplaceApproximation;
catch
RAW{4,icol}  = NaN;
end
offset=0;
for j=1:length(rmse),
    
  iraw = strmatch(deblank(lnam(j,:)),vnam,'exact');
  if ~isempty(iraw),  
      RAW{iraw+4,icol}=rmse(j);
  else
      offset=offset+1;
      RAW=[RAW(1:size(vnam,1)-1+offset+4-1,:);cell(1,size(RAW,2));RAW(size(vnam,1)-1+offset+4:end,:)];
      RAW{size(vnam,1)-1+offset+4,icol}=rmse(j);
      RAW{size(vnam,1)-1+offset+4,icol1+1}=deblank(lnam(j,:));
      RAW{size(vnam,1)-1+offset+4,icol+1}=deblank(lnam(j,:));
  end
end

offset1=offset;
for j=1:length(rmseK),
    
  iraw = strmatch(deblank(lnamK(j,:)),vnamK,'exact');
  if ~isempty(iraw),  
      RAW{iraw(1)+irawK+2+offset1,icol}=rmseK(j);
  else
      offset=offset+1;
      RAW=[RAW(1:size(vnam,1)+size(vnamK,1)-1+offset+5-1,:);cell(1,size(RAW,2));RAW(size(vnam,1)+size(vnamK,1)-1+offset+5:end,:)];
      RAW{size(vnam,1)+size(vnamK,1)-1+offset+5,icol}=rmse(j);
      RAW{size(vnam,1)+size(vnamK,1)-1+offset+5,icol1+1}=deblank(lnam(j,:));
      RAW{size(vnam,1)+size(vnamK,1)-1+offset+5,icol+1}=deblank(lnam(j,:));
  end
end

irawR2 = strmatch('r2',TXT(:,2));
RAW{irawR2+2+offset,icol}  = titcell;
offset2=offset;
for j=1:length(r2),
    
  iraw = strmatch(deblank(lnam(j,:)),vnam,'exact');
  
  if ~isempty(iraw),
      RAW{iraw+irawR2+2+offset2,icol}=r2(j);
  else
      offset=offset+1;
      RAW=[RAW(1:size(vnam,1)+size(vnamK,1)+size(vnam,1)-1+offset+6-1,:);cell(1,size(RAW,2));RAW(size(vnam,1)+size(vnamK,1)+size(vnam,1)-1+offset+6:end,:)];
      RAW{size(vnam,1)+size(vnamK,1)+size(vnam,1)-1+offset+6,icol}=rmse(j);
      RAW{size(vnam,1)+size(vnamK,1)+size(vnam,1)-1+offset+6,icol1+1}=deblank(lnam(j,:));
      RAW{size(vnam,1)+size(vnamK,1)+size(vnam,1)-1+offset+6,icol+1}=deblank(lnam(j,:));
  end
end
if ~ismac
    [SUCCESS,MESSAGE]=xlswrite(fname,RAW,sheet);
else
    [SUCCESS]=xlswrite_MACOS(fname,RAW,sheet);
end

