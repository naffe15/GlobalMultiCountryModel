function y0 = getIRFRuns(jxj,jexo,type)
global M_

if nargin<3, type='posterior'; end,

if strcmpi(type,'posterior')
  dirname = [CheckPath('metropolis') '/' ];
elseif strcmpi(type,'gsa')
  dirname = [CheckPath('GSA') '/' ];
else
  dirname = [CheckPath('prior') '/' ];
end  
filIRF = dir([dirname '/' M_.fname '_IRF_DSGEs*.mat']);

y0=[];
for file = 1:length(filIRF),
  load([dirname '/' M_.fname '_IRF_DSGEs' int2str(file)]);
%   y0 = [y0; squeeze(STOCK_IRF_DSGE(:,jxj,jexo,:))];
  y0 = [y0; (STOCK_IRF_DSGE(:,jxj,jexo,:))];
end
