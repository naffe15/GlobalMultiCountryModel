function m = max_nan(y);
% computes the mean of y with missing values

[nr,nc]=size(y);
if nr==1,
    y =y';
    nr=nc;
    nc = 1;
end

for j=1:nc,
    indx = find(~isnan(y(:,j)));
    m(1,j) = max(y(indx,j));
end
