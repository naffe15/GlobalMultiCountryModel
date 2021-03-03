%% Replication Package

% Alice Albonico, Ludivic Calés, Roberta Cardani, Olga Croitorov, Filippo Ferroni, 
% Massimo Giovannini, Stefan Hohberger, Beatrice Pataracchia, Filippo Pericoli, 
% Rafal Raciborski, Marco Ratto, Werner Roeger, Lukas Vogel

% Comparing post-crisis dynamics across Euro Area countries with the Global Multi-country model

% Economic Modelling 81(C), 2019, 242-273
%%
close all; clear all;

% add dynare path
addpath C:\dynare\4.5.6\matlab\

% Adding utilities
addpath util
addpath util_4.5b


% Specify output folders
workdirectory = pwd;
my_DE = 'DE\';
my_FR = 'FR\';
my_IT = 'IT\';
my_ES = 'ES\';
my_output = 'output_util\';

%% Specify country
DE = 1;
FR = 1;
IT = 1;
ES = 1;

%% Run DE model

cd(workdirectory)
if DE
    cd(my_DE)
    dynare gemc console nointeractive
end
cd(workdirectory)


%% Run FR model

cd(workdirectory)
if FR
    cd(my_FR)
    dynare gemc console nointeractive
end
cd(workdirectory)


%% Run IT model

cd(workdirectory)
if IT
    cd(my_IT)
    dynare gemc console nointeractive
end
cd(workdirectory)


%% Run ES model

cd(workdirectory)
if ES
    cd(my_ES)
    dynare gemc console nointeractive
end
cd(workdirectory)


%% Creating figures and plots in the paper

% Capacity utilization in the model and the data (Fig. 1)
cd(my_output)
CU_Figure1;

%%
% Dynamic responses (Figs. 2-7)
cd(my_output)
IRFs_Fig2_7;

%%
% Annual fit - 1- and 2-year ahead forecast (Figs. B.1 - B.4)
cd(my_output)
Fit_Fig_B1_B4;

%%
% Lead-lag structure of output growth and its main components (Figs. B.5 - B.8)
cd(my_output)
Lead_lag_Fig_B5_B8;

