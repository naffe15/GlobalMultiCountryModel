function name2cell(varargin);

cname = varargin(2:end);
if strcmp(varargin{1}(1),'-'),
  assignin('caller',varargin{1}(2:end),cname)
else
  disp('The name of the cell variable MUST start with - (minus)!!!')
end
