function data_moments = data_covariances(varargin);

% function data_covariances(varargin);
% compute moments in the data 
% ACF and ACF for the data are computed without  prefiltering, i.e.
% consitent with model estimation where data fluctuate around the  model
% steady state and not the sample mean
%
% inputs: list of endogenous
%
% outputs: 
% data_moments: wrap-up of data ACF/CCF for the declared list of
%     edogenous
%
global M_ oo_ options_

lambda=options_.hp_filter;
nlags = options_.ar;
fname_ = M_.fname;
[DirectoryName, info] = CheckPath('Output',fname_);
data_moments.var=zeros(length(varargin),length(varargin));
data_moments.corr=zeros(length(varargin),length(varargin));
data_moments.autocorr=zeros(length(varargin),nlags);
data_moments.acf=struct();
data_moments.ccf=struct();

[dataset_, dataset_info, newdatainterfaceflag] = makedataset(options_);
for j=1:length(varargin)
    data(:,j) = dataset_.data(:,j) - get_mean(varargin{j});
    if lambda,
        [s,desvabs] = hpfilter(data(:,j),lambda);
        data(:,j) = data(:,j) - s;
    end
    cova(j,j)=sum(data(:,j).^2)/size(data,1);
    std_dev(j) = sqrt(cova(j,j));
    cc(j,j)=1;
    for ii=1:j-1
        cova(j,ii)=sum(data(:,j).*data(:,ii))/size(data,1);
        cova(ii,j)=cova(j,ii);
        cc(j,ii)=cova(j,ii)/std_dev(j)/std_dev(ii);
        cc(ii,j) = cc(j,ii);
    end
end

data_moments.var=cova;
data_moments.corr=cc;
n=size(data,1);

for ii=-nlags:nlags,
    sigmaXCF(ii+nlags+1,1) = sqrt(1/(n-abs(ii)));
    BoundsXCF(ii+nlags+1,[1 2]) = sigmaXCF(ii+nlags+1,1)*[2 -2];
end
for j=1:length(varargin)
    % autocorr/crosscorr: 0: prefilter: assume data are without the mean,
    % consistent with moments computed around the steady state of  the
    % model
    if nlags,
        [ACF,Lags] = autocorr(data(:,j),nlags,0,2,0);
        for ii=1:nlags,
            sigmaQ(ii,1) = sqrt((1 + 2*(ACF(2:ii)'*ACF(2:ii)))/n);
            Bounds(ii,[1 2]) = sigmaQ(ii,1)*[2 -2];
        end
        data_moments.autocorr(j,:)=ACF(2:end);
        data_moments.acf = setfield(data_moments.acf,varargin{j},[Lags(2:end) ACF(2:end) Bounds]);
    end
    for i=1:length(varargin)
        [XCF,Lags] = crosscorr(data(:,j),data(:,i),nlags,2,0);
        data_moments.ccf = setfield(data_moments.ccf,[varargin{j} '_' varargin{i}],[Lags XCF BoundsXCF]);
    end
end



