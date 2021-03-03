function out=ForecastAnalysisFn(nfrcst,lastvalue,T, ccova, oo_ , M_, timebds)



% % % % % Input:
% % % % % nfrcst= number of forecasted quarters
% % % % % lastvalue: last quarter to include in forecast e.g.  2030.75
% % ss0 = getSmootherInfo(oo_);
%E_RPEDUM=0;
% % % % RPREME0=RPREME;
% % % % %RPREME=0.01;
% % % % % % ss_ = getSmootherInfo(oo_);
% % % % % % % % % % E_TAXDUM=0;
% % % % % % %RHORPEE0=RHORPEE;
% % % % % % %RHORPEE=0.99;
% % % % % % %TGOVB1E=0.0025;
% % % % % % %TGOVB1E=0.075;
% % % % % % %TGOVB2E=TGOVB1E*30;
% % % % % % steady;
% % % % % % check;
% % % % % % [Tnew,Rnew,SteadyState] = dynare_resolve;
% % % % % % yfcst=ss_.a;
% % % % % % for j=1:nfrcst,
% % % % % %     yfcst(:,j+1)=Tnew*yfcst(:,j);
% % % % % % end
% % % % % % yfcst0=ss_.a;
% % % % % % Vfrcst=zeros(length(Tnew), length(Tnew), nfrcst);
% % % % % % for j=1:nfrcst,
% % % % % %     yfcst0(:,j+1)=ss_.T*yfcst0(:,j);
% % % % % %     Vfrcst(:,:,j+1)=ss_.T*Vfrcst(:,:,j)*(ss_.T)'+ss_.R*ccova*(ss_.R)';
% % % % % % end
% % % % % % 
% % % % % % fnam = fieldnames(oo_.SmoothedVariables);
% % % % % % for j=1:length(fnam),
% % % % % %     eval(['ForecastedVariables.',fnam{j},'=[oo_.SmoothedVariables.',fnam{j},'(1:end); yfcst(j,2:end)''];'])
% % % % % %     eval(['ForecastedVariables0.',fnam{j},'=[oo_.SmoothedVariables.',fnam{j},'(1:end); yfcst0(j,2:end)''];'])
% % % % % %     eval(['ForecastedVariablesStd.',fnam{j},'=[oo_.SmoothedVariables.',fnam{j},'(1:end)*NaN; squeeze(sqrt(Vfrcst(j,j,2:end)))];'])
% % % % % % end

outFV=ForecastedVariablesFn(nfrcst,ccova,oo_);
ForecastedVariables=outFV.ForecastedVariables;
% ForecastedVariables0=outFV.ForecastedVariables0;
ForecastedVariablesStd=outFV.ForecastedVariablesStd;
% % % % % % % % /*
% % % % % % % % pplotvar= {'E_GY', 'E_GIHOUSE', 'E_GIHOUSECC', 'E_INFHOUSE', 'E_INFLAND', 'E_GC', 'E_GCCC', 'E_INFC', 'E_INFY', 'E_GL', 'E_GIG', 'E_GI','E_R',...
% % % % % % % % 'E_BY', 'E_L','E_INOM','E_GEX','E_GIM','E_GK','E_GBY', 'E_GBY_Q', 'E_WRINF', 'E_TBYN','E_LCY','E_LIY','E_LIHOUSEY','E_LICONSTRY',...
% % % % % % % % 'E_LHOUSEY','E_LKY','E_LCSN','E_LISN','E_LIHOUSESN','E_LICONSTRSN','E_LHOUSESN','E_GG', 'E_TOT'};
% % % % % % % % */



% % % % % % % /*pplotvar= {'E_GY',...
% % % % % % % 'E_GULC', ...
% % % % % % % 'E_GTFP', ...
% % % % % % % 'E_GUPTOT', ...
% % % % % % % 'E_GBY',...
% % % % % % % 'E_PHI', ...
% % % % % % % 'E_LCY',...
% % % % % % % 'E_LIY',...
% % % % % % % 'E_LICONSTRY',...
% % % % % % % 'E_L',...
% % % % % % % 'E_TBYN',...
% % % % % % % 'E_INFY', ...
% % % % % % % 'E_REER',...
% % % % % % % 'E_LGY', ...
% % % % % % % 'E_LIGY', ...
% % % % % % % 'E_LTRY',...
% % % % % % % 'E_TL',...
% % % % % % % 'E_BY', ...
% % % % % % % 'E_LTFP', ...
% % % % % % % 'E_LUPTOT'};
% % % % % % % */

% % % % % % % % /*1.
% % % % % % % % GY    GTFP    GUPTOT
% % % % % % % % LCY   LIY       LICONSTRY
% % % % % % % % L       PHI
% % % % % % % %
% % % % % % % % 2.
% % % % % % % % GBY (Quarterly!)   BY     TL
% % % % % % % % LGY                   LIGY     LTRY
% % % % % % % % TBYN             REER       GULC
% % % % % % % % */
pplotvar= {'E_GBY_Q',...
    'E_BY', ...
    'E_TL',...
    'E_LGY', ...
    'E_LIGY', ...
    'E_LTRY',...
    'E_TBYN',...
    'E_REER',...
    'E_GULC',...
    'E_GY',...
    'E_GTFP', ...
    'E_GUPTOT', ...
    'E_LCY',...
    'E_LIY',...
    'E_LICONSTRY',...
    'E_L',...
    'E_PHI'};



GPOPpos=strmatch('W_GPOP', M_.param_names);
GP0pos=strmatch('GP0', M_.param_names);
GTFPpos=strmatch('W_GTFP', M_.param_names);

W_GPOP=M_.params(GPOPpos);
W_GTFP=M_.params(GTFPpos);
GP0=M_.params(GP0pos);

ForecastedVariables1=ForecastedVariables;
ForecastedVariables1.E_GY=ForecastedVariables1.E_GY+W_GPOP;
ForecastedVariables1.E_GULC=ForecastedVariables1.E_GULC+GP0-W_GTFP-W_GPOP;
ForecastedVariables1.E_GTFP=ForecastedVariables1.E_GTFP+W_GTFP+W_GPOP;
ForecastedVariables1.E_GUPTOT=ForecastedVariables1.E_GUPTOT+W_GTFP+W_GPOP;



SmoothedVariables0=oo_.SmoothedVariables;
SmoothedVariables0.E_GY=SmoothedVariables0.E_GY+W_GPOP;
SmoothedVariables0.E_GULC=SmoothedVariables0.E_GULC+GP0-W_GTFP-W_GPOP;
SmoothedVariables0.E_GTFP=SmoothedVariables0.E_GTFP+W_GTFP+W_GPOP;
SmoothedVariables0.E_GUPTOT=SmoothedVariables0.E_GUPTOT+W_GTFP+W_GPOP;

for j=1:length(pplotvar)
    steadyvar(j)=get_mean(pplotvar{j});
end

steadyvar(10)=steadyvar(10)+W_GPOP;
steadyvar(9)=steadyvar(9)+GP0-W_GTFP-W_GPOP;
steadyvar(11)=steadyvar(11)+W_GTFP+W_GPOP;
steadyvar(12)=steadyvar(12)+W_GTFP+W_GPOP;



TT=1995:0.25:lastvalue;
jfig=0;
for j=1:length(pplotvar)
    if mod(j,9)==1
        figure, jplot=0; jfig=jfig+1;
    end
    jplot=jplot+1;
    subplot(3,3,jplot)
    %vplot=eval(['ForecastedVariables11.',pplotvar{j},';']);
    vplot0=eval(['ForecastedVariables1.',pplotvar{j},';']);
    vplotStd1=eval(['ForecastedVariables1.',pplotvar{j},'+2*ForecastedVariablesStd.',pplotvar{j},';']);
    vplotStd_1=eval(['ForecastedVariables1.',pplotvar{j},'-2*ForecastedVariablesStd.',pplotvar{j},';']);
    % % % % % % % % % % % % /*vplot0=eval(['ForecastedVariables.',pplotvar{j},';']);
    % % % % % % % % % % % % vplotStd1=eval(['ForecastedVariables.',pplotvar{j},'+2*ForecastedVariablesStd.',pplotvar{j},';']);
    % % % % % % % % % % % % vplotStd_1=eval(['ForecastedVariables.',pplotvar{j},'-2*ForecastedVariablesStd.',pplotvar{j},';']);
    % % % % % % % % % % % % */
    
    hold on, plot(TT,vplot0(1:length(TT))+get_mean(pplotvar{j}),'b'),
    hold on, plot(TT,[vplotStd1(1:length(TT)) vplotStd_1(1:length(TT))]+get_mean(pplotvar{j}), 'b:')
    title(pplotvar{j},'interpreter','none')
    %vplot=eval(['oo_.SmoothedVariables.',pplotvar{j},';']);
    vplot=eval(['SmoothedVariables0.',pplotvar{j},';']);
    hold on, plot(T,vplot+get_mean(pplotvar{j}),'r.')
    %hold on, plot([TT(1) TT(end)],[get_mean(pplotvar{j}) get_mean(pplotvar{j})],'k:')
    
    %hold on, plot([TT(1) TT(end)],[get_mean(pplotvar{j}) get_mean(pplotvar{j})],'k:')
    hold on, plot([TT(1) TT(end)],[steadyvar(j) steadyvar(j)],'k:')
if isempty(timebds) || nargin<7
    set(gca,'xlim',[TT(1)-1 TT(end)+0.25]);
else
    set(gca,'xlim',[timebds(1) timebds(2)]);
end
    if jplot==9 || j==length(pplotvar)
        saveas(gcf,[M_.fname,'_smooth_and_frcst_short',int2str(jfig)]);
        %saveas(gcf,[M_.fname,'_smooth_and_frcst_',int2str(jfig)]);
    end
end

% out.ForecastedVariables0=ForecastedVariables0;
out.ForecastedVariables=ForecastedVariables;

return

