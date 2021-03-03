function [var, ny, nx, forecast_data] = vecm_toolbox(nlags)
% vecm_toolbox  Routines for VECM methods
%
% [var, ny, nx, forecast_data] = vecm_toolbox(nlags)
%
% Computes several things for the estimations of a VECM(nlags):
% ny:         number of endogenous variables
% nx:         number of exogenous variables (equal to zero, or one if a
%             constant term is included)
%
% This function uses the following Dynare options:
% - datafile, first_obs, varobs, xls_sheet, xls_range, nobs, presample, nk
% - opts_vecm_{varobs,vecm}
    
    global options_ opts_vecm_
    
    % Load dataset
    dataset = read_variables(options_.datafile, opts_vecm_.varobs, [], [], []);
    datax = read_variables(options_.datafile, opts_vecm_.vecm, [], [], []);
    try
      TT = read_variables(options_.datafile, 'T', [], [], []);
    catch
      TT=[1:size(dataset,1)]';
    end
    first_obs = max(2,options_.first_obs);
%     options_ = set_default_option(options_, 'nobs', size(dataset,1)-first_obs);
    nobs = size(dataset,1)-first_obs+1;
    % Parameters for vecm
    nk=options_.nk;
    for j=1:size(datax,2),
      yyy=[diff(datax(:,j))];
      beta(j,:) = (inv(dataset(first_obs:end,:)'*dataset(first_obs:end,:))*dataset(first_obs:end,:)'*yyy(first_obs-1:end,:))';
    end

    
    if first_obs + options_.presample <= nlags
        error('first_obs+presample should be > nlags (for initializing the VAR)')
    end

    idx = first_obs+options_.presample-nlags:first_obs+nobs-1;
    
    % Prepare dataset
    if options_.loglinear & ~options_.logdata
        dataset = log(dataset);
        datax = log(datax);
    end
     
    ny = size(dataset, 2);
    nx = size(datax,2);
    if options_.prefilter | options_.noconstant,
        % use sample mean to de-mean the data
        y0=mean(dataset(idx,:));
        x0=mean(datax(idx-1,:));
        dataset = dataset - ones(length(TT),1)*y0;
        datax = datax - ones(length(TT),1)*x0;
    else
        % use model steady state to de-mean the data
      for j=1:size(dataset,2),
        y0(j)=get_mean(opts_vecm_.varobs(j,:));
        dataset(:,j)=dataset(:,j)-y0(j);        
      end
      for j=1:size(datax,2),
        x0(j)=get_mean(opts_vecm_.vecm(j,:));
        datax(:,j)=datax(:,j)-x0(j);        
      end
        
    end
    

    ydata = dataset(idx, :);
    T = size(ydata, 1);
    if (size(ydata,2)*nlags+size(datax,2))>(T-nlags),
      warning('Not enough degree of freedhom for VECM, use BVAR instead!')
      var=[];
      if nargout >= 4
        forecast_data=[];    
      end
      return
    end
    xdata = datax(idx-1, :);

    % VECM estimation
    var = rfvar3(ydata, nlags, xdata, round(beta), nk);
    var.ydata=ydata;
    var.T=TT(idx);
    Tu = size(var.u, 1);
    
    var.S = var.u(:,:,1)' * var.u(:,:,1);
    var.beta=round(beta);
    
    
    % Add forecast informations
    if nargout >= 4
        forecast_data.initval = ydata(end-nlags+1:end, :);
        forecast_data.xval = xdata(end, :);
        forecast_data.x0=x0;
        forecast_data.y0=y0;
        if first_obs + nobs <= size(dataset, 1)
            forecast_data.realized_val = dataset(first_obs+nobs:end, :);
            forecast_data.realized_xdata = datax(first_obs+nobs:end, :);
        else
            forecast_data.realized_val = [];
        end
    end
    
  

function var=rfvar3(ydata,lags,xdata,beta,nk,breaks)
%function var=rfvar3(ydata,lags,xdata,beta,nk,breaks)
% This algorithm goes for accuracy without worrying about memory requirements.
% ydata:   dependent variable data matrix
% xdata:   exogenous variable data matrix
% lags:    number of lags
% breaks:  rows in ydata and xdata after which there is a break.  This allows for
%          discontinuities in the data (e.g. war years) and for the possibility of
%          adding dummy observations to implement a prior.  This must be a column vector.
%          Note that a single dummy observation becomes lags+1 rows of the data matrix,
%          with a break separating it from the rest of the data.  The function treats the 
%          first lags observations at the top and after each "break" in ydata and xdata as
%          initial conditions. 
%      The program assumes that the first lags rows of ydata and xdata are real data, not dummies.
%      Dummy observations should go at the end, if any.  If pre-sample x's are not available,
%      repeating the initial xdata(lags+1,:) row or copying xdata(lags+1:2*lags,:) into 
%      xdata(1:lags,:) are reasonable subsititutes.  These values are used in forming the
%      persistence priors.
% Adapted from code written by Christopher Sims.  This version 6/15/03.
    [T,nvar] = size(ydata);
    nox = isempty(xdata);
    if ~nox
        [T2,nx] = size(xdata);
    else
        T2 = T;
        nx = 0;
        xdata = zeros(T2,0);
    end
    % note that x must be same length as y, even though first part of x will not be used.
    % This is so that the lags parameter can be changed without reshaping the xdata matrix.
    if T2 ~= T, error('Mismatch of x and y data lengths'),end
    if nargin < 5,
      nk=4;
    end
    if nargin < 6
        nbreaks = 0;
        breaks = [];
    else
        nbreaks = length(breaks);
    end
    breaks = [0;breaks;T];
    smpl = [];
    for nb = 1:nbreaks+1
        smpl = [smpl;[breaks(nb)+lags+1:breaks(nb+1)]'];
    end
    Tsmpl = size(smpl,1);
    X = zeros(Tsmpl,nvar,lags);
    for is = 1:length(smpl)
        X(is,:,:) = ydata(smpl(is)-(1:lags),:)';
    end
    X = [X(:,:) xdata(smpl,:)];
    y = ydata(smpl,:);
    % Everything now set up with input data for y=Xb+e 
    
    
    % Compute OLS regression and residuals
    [vl,d,vr] = svd(X,0);
    di = 1./diag(d);
    B = (vr.*repmat(di',nvar*lags+nx,1))*vl'*y;
    u(:,:,1) = y-X*B;
    xxi = vr.*repmat(di',nvar*lags+nx,1);
    xxi = xxi*xxi';
    
    yhat(:,:,1)=X*B;
    x0=xdata(smpl,:);
    xhat=x0+yhat(:,:,1)*beta';
    for j=2:nk,
      X=[yhat(:,:,j-1) X(:,1:nvar*(lags-1)) xhat];
      yhat(:,:,j)=X*B;
      xhat=xhat+yhat(:,:,j)*beta';
      u(j:end,:,j)=y(j:end,:)-yhat(1:end-j+1,:,j);
    end
    
    var.B = B;
    var.u = u;
    var.xxi = xxi;
    var.yhat = yhat;
