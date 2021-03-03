# GlobalMultiCountryModel
Replication Package
***************************************************************************

Economic Modelling 81(C), 2019, 242-273
Title: Comparing post-crisis dynamics across Euro Area countries with the Global Multi-country model
Authors: Alice Albonico, Ludivic Calés, Roberta Cardani, Olga Croitorov, Filippo Ferroni, 
Massimo Giovannini, Stefan Hohberger, Beatrice Pataracchia, Filippo Pericoli, 
Rafal Raciborski, Marco Ratto, Werner Roeger, Lukas Vogel

Corresponding author: Marco Ratto

***************************************************************************

The package has been tested with:
MATLAB 2017b
Dynare 4.5.6 (https://www.dynare.org/release/) 

Note: This package should be compatible with other Dynare 4.5 versions, 
but may not run properly with newer Dynare versions (4.6 onward)
  
***************************************************************************


To replicate the paper results, please proceed as follows:

1)  Add your Dynare path to the Matlab search path:

        addpath [YOUR DYNARE PATH] \dynare\4.5.6\matlab\

2)  Main_replic.m runs and replicates all results for all four countries:
    It adds the utilities, loads the data and the posterior mode (gemc_auto_mode) [i.e. Table 2] 
    and generates the results for all four country versions:
    - DE
    - FR
    - IT
    - ES 
    Note: The default Main_replic.m runs ALL four countries sequentially.

3)  All results are saved in the respective model subfolders:
    - gemc_moments.xls (Theoretical moments and model fit, see Table 3)
    - gemc_annual_fit_1.fig (Unconditional 1- and 2-year ahead forecast, see Figures B.1-B.4)
    - gemc\graphs (shock decompositions of yoy real GDP growth and the qoq trade balance-to-GDP ratio, see Figures 8-15)
    - gemc\Output (Lead-lag structure of output growth and its main components, see Figures B.5-B.8)
    - gemc\Output_paper (Impulse response functions, see Figures 2-7)

4)  The main figures can be also reproduced by the routines labelled as in the paper (e.g. CU_Figure1.m).
    Note: You have to run the model simulations first in order to have all the results saved.


All country versions are nested and can be run individually by setting the country options to 0 or 1. 
For example, if you want to run only the Italian model, please specify the model options as:
- DE = 0 
- FR = 0
- IT = 1
- ES = 0

Additionally, if you want to run a model without the Main_replic.m script, please proceed as follows:
- Add the utilities and Dynare to the path
    addpath util
    addpath util_4.5
    addpath [YOUR DYNARE PATH] \dynare\4.5.6\matlab 
- Go into the respective country subfolder (DE, FR, IT, ES)
- Run the model by typing in the command window:
    dynare gemc console nointeractive

