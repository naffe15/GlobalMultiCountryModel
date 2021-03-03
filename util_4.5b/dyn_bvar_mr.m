function dyn_bvar(maxnlags,tau,d,lambda,mu,omega,train,flat,breaks)
% stephane.adjemian@cepremap.cnrs.fr [21 april, 2004]


global options_

eval(options_.datafile);
dataset = [ ];

for i=1:size(options_.varobs,1)
    dataset = [dataset eval(deblank(options_.varobs(i,:)))];
end    

default_nlags  = 8;
default_const  = 1;
default_tau    = 3;
default_d      = 0.5;
default_lambda = 5;
default_mu     = 2;
default_omega  = 1;

if nargin == 0
    maxnlags =  default_nlags;
    tau      =  default_tau;
    d        =  default_d;
    lambda   =  default_lambda;    
    mu       =  default_mu;
    omega    =  default_omega;
    breaks   = [];
    train    = [];
    flat     = 0;
elseif nargin == 1
    tau    =  default_tau;
    d      =  default_d;
    lambda =  default_lambda;    
    mu     =  default_mu;
    omega  =  default_omega;
    breaks = [];
    train  = [];
    flat   = 0;
elseif nargin == 2
    if isempty(maxnlags)
        maxnlags = default_nlags;
    end   
    d      =  default_d;
    lambda =  default_lambda;    
    mu     =  default_mu;
    omega  =  default_omega;
    breaks = [];
    train  = [];
    flat   = 0;
elseif nargin == 3
    if isempty(maxnlags)
        maxnlags = default_nlags;
    end   
    if isempty(tau)
        tau = default_tau;
    end   
    lambda =  default_lambda;    
    mu     =  default_mu;
    omega  =  default_omega;
    breaks = [];
    train  = [];
    flat   = 0;
elseif nargin == 4
    if isempty(maxnlags)
        maxnlags = default_nlags;
    end   
    if isempty(tau)
        tau   = default_tau;
    end
    if isempty(d)
        d     = default_d;
    end   
    mu     =  default_mu;
    omega  =  default_omega;
    breaks = [];
    train  = [];
    flat   = 0;
elseif nargin == 5
    if isempty(maxnlags)
        maxnlags = default_nlags;
    end   
    if isempty(tau)
        tau   = default_tau;
    end   
    if isempty(d)
        d     = default_d;
    end   
    if isempty(lambda)
        lambda =  default_lambda;    
    end    
    omega  =  default_omega;
    breaks = [];
    train  = [];
    flat   = 0;
elseif nargin == 6
    if isempty(maxnlags)
        maxnlags = default_nlags;
    end   
    if isempty(tau)
        tau = default_tau;
    end   
    if isempty(d)
        d = default_d;
    end   
    if isempty(lambda)
        lambda = default_lambda;
    end   
    if isempty(mu)
        mu = default_mu;
    end   
    train  = [];
    flat   = 0;
    breaks = [];
elseif nargin == 7
    if isempty(maxnlags)
        maxnlags = default_nlags;
    end   
    if isempty(tau)
        tau = default_tau;
    end   
    if isempty(d)
        d = default_d;
    end   
    if isempty(lambda)
        lambda = default_lambda;
    end   
    if isempty(mu)
        mu = default_mu;
    end
    if isempty(omega)
        omega = default_omega;
    end   
    flat   = 0;
    breaks = [];
elseif nargin == 8
    if isempty(maxnlags)
        maxnlags = default_nlags;
    end   
    if isempty(const)
        const = default_const;
    end   
    if isempty(tau)
        tau = default_tau;
    end   
    if isempty(d)
        d = default_d;
    end   
    if isempty(lambda)
        lambda = default_lambda;
    end   
    if isempty(mu)
        mu = default_mu;
    end
    if isempty(omega)
        omega = default_omega;
    end
    breaks = [];
    train  = [];
elseif nargin == 9
    if isempty(maxnlags)
        maxnlags = default_nlags;
    end   
    if isempty(const)
        const = default_const;
    end   
    if isempty(tau)
        tau = default_tau;
    end   
    if isempty(d)
        d = default_d;
    end   
    if isempty(lambda)
        lambda = default_lambda;
    end   
    if isempty(mu)
        mu = default_mu;
    end
    if isempty(omega)
        omega = default_omega;
    end
    train = [];
    flat  = 0;
elseif nargin > 9
    disp('dyn_bvar :: too many arguments.')
end

%load cul;
%dataset=xx;
% tried with VECM!!!!

mnprior.tight = tau;
mnprior.decay = d;
vprior.sig = std(dataset(options_.first_obs+options_.presample-maxnlags:options_.first_obs+options_.presample-1,:))';
vprior.w = omega;
xdata=[];

for lag = 1:maxnlags
    ydata = dataset(options_.first_obs+options_.presample-lag:options_.first_obs+options_.nobs-1,:);
    xdata = xvec(options_.first_obs+options_.presample-lag:options_.first_obs+options_.presample+options_.nobs-1,:);
    w=mgnldnsty(ydata,lag,[ones(options_.nobs+lag,1) xdata],breaks,lambda,mu,mnprior,vprior,train,flat);
    disp(' ')
    fprintf('The marginal log density of the BVAR(%g) model is equal to %10.4f \n',lag,w);
    disp(' ')
end

