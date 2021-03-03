function [b, bstd, sig, res] = regre_nan(x,y)
% function [b, bstd, sig, res] = regre_nan(x,y),
% performs a regression neglecting missing values

if nargin == 0,
    disp('[b, bstd, sig, res] = regre_nan(x,y);')
    return
end
y=y(:);
if size(x,1)~=length(y),
    error('x dimension mismatch w.r.t. y');
end

indy = find(~isnan(sum([x y],2)));
b=x(indy,:)\y(indy);
res = y(indy)-x(indy,:)*b;
sig = sqrt(res'*res/(length(indy)-length(b)));
bstd = sig*sqrt(diag(inv(x(indy,:)'*x(indy,:))));
