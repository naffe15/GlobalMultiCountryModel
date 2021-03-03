function  plot_forecasts(pplotvar, bayestopt_, M_,options_, oo_, ForecastedVariables0,ForecastedVariablesStd,GP0, W_GPOP, W_GTFP,T )



% function mode_check(x,fval,hessian,gend,data,lb,ub)
% Checks the maximum likelihood mode
%
% INPUTS
%    x:       mode
%    fval:    value at the maximum likelihood mode
%    hessian: matrix of second order partial derivatives
%    gend:    scalar specifying the number of observations
%    data:    matrix of data
%    lb:      lower bound
%    ub:      upper bound
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2010 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.





ForecastedVariables1=ForecastedVariables0;
ForecastedVariables1.E_GY=ForecastedVariables1.E_GY+W_GPOP;
ForecastedVariables1.E_GTFP=ForecastedVariables1.E_GTFP+W_GTFP+W_GPOP;
ForecastedVariables1.E_GULC=ForecastedVariables1.E_GULC+GP0-W_GTFP-W_GPOP;
ForecastedVariables1.E_GUPTOT=ForecastedVariables1.E_GUPTOT+W_GTFP+W_GPOP;



SmoothedVariables0=oo_.SmoothedVariables;
SmoothedVariables0.E_GY=SmoothedVariables0.E_GY+W_GPOP;
SmoothedVariables0.E_GTFP=SmoothedVariables0.E_GTFP+W_GTFP+W_GPOP;
SmoothedVariables0.E_GULC=SmoothedVariables0.E_GULC+GP0-W_GTFP-W_GPOP;
SmoothedVariables0.E_GUPTOT=SmoothedVariables0.E_GUPTOT+W_GTFP+W_GPOP;

for j=1:length(pplotvar)
    steadyvar(j)=get_mean(pplotvar{j});
    if strcmp(pplotvar{j}, 'E_GY'),
        steadyvar(j)=steadyvar(j)+W_GPOP;
    end
    if strcmp(pplotvar{j}, 'E_GTFP'),
        steadyvar(j)=steadyvar(j)+W_GTFP+W_GPOP;
    end
    if strcmp(pplotvar{j}, 'E_GUPTOT'),
        steadyvar(j)=steadyvar(j)+W_GTFP+W_GPOP;
    end
    if strcmp(pplotvar{j}, 'E_GULC'),
        steadyvar(j)=steadyvar(j)+GP0-W_GTFP-W_GPOP;
    end
        
end

% steadyvar(10)=steadyvar(10)+W_GPOP;
% steadyvar(11)=steadyvar(11)+W_GTFP+W_GPOP;
% steadyvar(9)=steadyvar(9)+GP0-W_GTFP-W_GPOP;
% steadyvar(12)=steadyvar(12)+W_GTFP+W_GPOP;


TT=1995:0.25:2030.75;
jfig=0;
for j=1:length(pplotvar)
    
    varnamej=pplotvar{j};
    varpos= strmatch(varnamej, M_.endo_names, 'exact');
    texname{j,:}=M_.endo_names_tex(varpos,:);
    varname{j,:}=varnamej;
    
    
    if mod(j,9)==1
        figure, jplot=0; jfig=jfig+1;
    end
    jplot=jplot+1;
    subplot(3,3,jplot)
    % //vplot=eval(['ForecastedVariables11.',pplotvar{j},';']);
    vplot0=eval(['ForecastedVariables1.',pplotvar{j},';']);
    vplotStd1=eval(['ForecastedVariables1.',pplotvar{j},'+2*ForecastedVariablesStd.',pplotvar{j},';']);
    vplotStd_1=eval(['ForecastedVariables1.',pplotvar{j},'-2*ForecastedVariablesStd.',pplotvar{j},';']);
    % /*vplot0=eval(['ForecastedVariables0.',pplotvar{j},';']);
    % vplotStd1=eval(['ForecastedVariables0.',pplotvar{j},'+2*ForecastedVariablesStd.',pplotvar{j},';']);
    % vplotStd_1=eval(['ForecastedVariables0.',pplotvar{j},'-2*ForecastedVariablesStd.',pplotvar{j},';']);
    % */
    
    hold on, plot(TT,vplot0(1:length(TT))+get_mean(pplotvar{j}),'b'),
    hold on, plot(TT,[vplotStd1(1:length(TT)) vplotStd_1(1:length(TT))]+get_mean(pplotvar{j}), 'b:')
    title(pplotvar{j},'interpreter','none')
    % //vplot=eval(['oo_.SmoothedVariables.',pplotvar{j},';']);
    vplot=eval(['SmoothedVariables0.',pplotvar{j},';']);
    hold on, plot(T,vplot+get_mean(pplotvar{j}),'r.')
    % //hold on, plot([TT(1) TT(end)],[get_mean(pplotvar{j}) get_mean(pplotvar{j})],'k:')
    
    % //hold on, plot([TT(1) TT(end)],[get_mean(pplotvar{j}) get_mean(pplotvar{j})],'k:')
    hold on, plot([TT(1) TT(end)],[steadyvar(j) steadyvar(j)],'k:')
    set(gca,'xlim',[TT(1)-1 TT(end)+0.25])
    if jplot==9 || j==length(pplotvar)
        saveas(gcf,[M_.fname,'_smooth_and_frcst_',int2str(jfig)]);
        eval(['print -dpdf ',M_.fname,'_smooth_and_frcst_',int2str(jfig),'.pdf'])
        eval(['print -depsc2 ',M_.fname,'_smooth_and_frcst_',int2str(jfig),'.eps'])
        
    end
end



fidTeX = fopen([M_.fname '_Forecasts.TeX'],'w');


files=dir('*frcst*.eps');
h1=1;
    
numberofvars=length(pplotvar);
for jj = 1:length(files),
    
    
    frcstfilename=files(jj).name;
    [a,b,c]=fileparts(frcstfilename);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    for kk=1:min(numberofvars,9)
%         hh=9*(h1-1)+1;
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(varname{h1,:}),['$' deblank(texname{h1,:}) '$']);
        h1=h1+1;
    end
    
    numberofvars=numberofvars-9;
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,['\\includegraphics[width=0.80\textwidth] {' b '} \n']);
    fprintf(fidTeX,'\\caption{Forecast}');
    fprintf(fidTeX,'\\label{Fig:Forecast:%s}\n',int2str(jj));
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,' \n');
end

fclose(fidTeX);


