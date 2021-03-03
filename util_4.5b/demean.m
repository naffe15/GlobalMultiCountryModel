function x=demean(x)

mx=mean_nan(x);
for j=1:size(x,2),
    x(:,j)=x(:,j)-mx(j);
end
