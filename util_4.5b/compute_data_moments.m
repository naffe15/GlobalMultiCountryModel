function [data_std, data_corr, data_autocorr] = compute_data_moments(oo, ar,varargin)

% persistent h;
T = length(getfield(oo.SmoothedVariables,varargin{1}));
for j=1:length(varargin)
 y(:,j) = getfield(oo.SmoothedVariables,varargin{j});
 data_var(j,j) = (sum(y(:,j).^2)./(T-1));
 data_std(j) = sqrt(data_var(j,j));
 data_corr(j,j)=1;
 for i=1:(j-1)
    data_corr(i,j) = sum(y(:,i).*y(:,j))/(T-1)./(data_std(j)*data_std(i));
    data_corr(j,i) = data_corr(i,j);
 end
 for k=1:ar
     ylagged(:,j)=lagged(y(:,j),k);
     data_autocorr{k}(j,j) = sum(ylagged(k+1:end,j).*y(k+1:end,j))/(T-k-1)./(data_std(j)*data_std(j));
 for i=1:(j-1)
    data_autocorr{k}(i,j) = sum(ylagged(k+1:end,j).*y(k+1:end,i))/(T-k-1)./(data_std(i)*data_std(j));
    data_autocorr{k}(j,i) = sum(ylagged(k+1:end,i).*y(k+1:end,j))/(T-k-1)./(data_std(i)*data_std(j));
 end
 end
end