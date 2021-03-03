function make_bvardata_mat(varargin)
global oo_ M_ options_

% This function generates a matfle contaitng the times series used for the
% BVAR estimation and forecasts. The time series are the variables obtained
% from smoothing the endogenous variables of the DSGE state space. 
% date:     10/02/2016
% author:   FF

string = ['save dataobs_bvar '];
for jj = 1 : length(varargin{1})
    temp = get_smooth(varargin{1}{jj});
    eval([varargin{1}{jj} '= temp(1:options_.nobs);'  ])    
    string = [ string ' ' varargin{1}{jj}];
end

eval(string);
