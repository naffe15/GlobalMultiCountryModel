function plot_Unobserved_Variables_Fn(T,oo_, options_,M_,modeltype)


if nargin<5 || isempty(modeltype)
    modeltype=1;
end

figure('name','Unobserved variables 1')
subplot(3,3,1),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_LTFP);
title('Labour Augmenting Technology Shock')
subplot(3,3,2),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_LUPI);
title('Investment Productivity Shock')
subplot(3,3,3)
  plot(T(1:options_.nobs),oo_.SmoothedShocks.E_EPS_M);
title('Monetary Shock')
subplot(3,3,4),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_RPREMK);
title('Shock to Stock Market Risk Premium')
subplot(3,3,5),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_RPREME);
title('Shock to International Capital Flows')
subplot(3,3,6),
  plot(T(1:options_.nobs),oo_.SmoothedShocks.E_EPS_XW+oo_.SmoothedShocks.E_EPS_PW);
title('Shock to External Demand')
subplot(3,3,7),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_L);
title('Shock to Wages')

subplot(3,3,8),
if modeltype
 plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_DEBTCCT);
title('Shock to Lending Conditions')
end
subplot(3,3,9),
if modeltype
    plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_RPREMHOUSENLC);
title('Shock to House Price Risk Premium')
end


saveas(gcf,'quest3hlmr_Unobserved_Variables_1.fig')
print -dpdf quest3hlmr_Unobserved_Variables_1.pdf
print -depsc2 quest3hlmr_Unobserved_Variables_1.eps


figure('name','Unobserved variables 2')
subplot(3,3,1),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_ETA);
title('Shock to markup domestic prices')
subplot(3,3,2),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_ETAM);
title('Shock to markup import prices')
subplot(3,3,3)
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_ETAX);
title('Shock to markup export prices')
subplot(3,3,4),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_G);
title('Shock to Goverment Consumption')
subplot(3,3,5),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_IG);
title('Shock to Goverment Investment')
subplot(3,3,6),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_TR);
title('Shock to Goverment Transfer')
subplot(3,3,7),
  plot(T(1:options_.nobs),oo_.SmoothedVariables.E_ZEPS_CNLC);
title('Shock to Confidence')
saveas(gcf,'quest3hlmr_Unobserved_Variables_2.fig')
print -dpdf quest3hlmr_Unobserved_Variables_2.pdf
print -depsc2 quest3hlmr_Unobserved_Variables_2.eps


%     if isempty(NAMES)
%         NAMES = name;
%         TeXNAMES = texname;
%     else
%         NAMES = char(NAMES,name);
%         TeXNAMES = char(TeXNAMES,texname);
%     end
%     
    fidTeX = fopen([M_.fname '_Unobserved_Variables.TeX'],'w+');


files=dir('*_Unobserved_Variables*.eps');

for jj = 1:length(files),
    
    unobservedfilename=files(jj).name;
    [a,b,c]=fileparts(unobservedfilename);
%     quest3hlmr_shock_dec_qoqt_E_LCY_Init
%     pos=findstr(b,'E_');
    
        % TeX eps loader file
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        %             fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,['\\includegraphics[width=0.80\textwidth] {' b '} \n']);
        fprintf(fidTeX,'\\caption{Unobserved Variables}');
        fprintf(fidTeX,'\\label{Fig: Unobserved Variables:%s}\n',int2str(jj));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
end


fclose(fidTeX);
