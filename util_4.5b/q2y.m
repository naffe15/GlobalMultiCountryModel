function y_a = q2y(y,flow)
% transforms quarterly data to yearly: 
% y_a = q2y(y,flow)
% flow = 1 (default)
% flow = 0 means stock variables

y=y(:);
if nargin<2,
    flow=1;
end
y_a = [NaN(3,1); y(1:end-3)+y(2:end-2)+y(3:end-1)+y(4:end)];
if flow==0,
    y_a=y_a./4;
end

