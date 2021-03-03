function check_mean(varargin)
global M_ oo_ options_

if exist(options_.datafile)
  instr = options_.datafile;
else
  instr = ['load ' options_.datafile];
end
eval(instr);
if ~exist('T','var'), 
  temp = eval(deblank(options_.varobs(1,:)));
  T=[1:length(temp)]'; 
  clear temp;
end
presample = options_.presample;

nobs = length(varargin);
fobs = options_.first_obs;

ifig=0;
iplo=0;
for j=1:nobs,
    iplo=iplo+1;
    if iplo>9,
        iplo=1;
    end
    if iplo==1,
        figure,
        ifig=ifig+1;
        
    end
    subplot(3,3,iplo),
%    try
%        vv=eval(varargin{j});
%        mm=get_mean(varargin{j});
%    catch
       vv=eval(['oo_.UpdatedVariables.',varargin{j}]);
       mm=get_mean(varargin{j});
%    end
        
    plot(T(fobs+presample:end),vv(fobs+presample:end),'k')
    hold on, plot(T([fobs,end]),[mm mm],'r'),
    set(gca,'xlim',T([fobs,end])), title(varargin{j},'interpreter','none')
    
end
