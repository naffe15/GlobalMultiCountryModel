function c = corrcoef_nan(y);
% computes the cov of y with missing values

[nr,nc]=size(y);
if nr==1,
    y =y';
    nr=nc;
    nc = 1;
end

indx = find(~isnan(sum(y,2)));
c = corrcoef(y(indx,:));
