function re_scaling(oo_, xname, xscale, vargin)

%E_Y E_C E_I E_CNLC E_CLC E_TBYN E_G E_IG E_TR E_L E_WR E_DBGYN E_LYGAP E_PHI E_PHIC;
% TR scaling
% WS = 0.4860
% TRW = 0.4023
% TRYN = 0.1955


nfig=0;
nplo=0;
for j=1:length(vargin),

  name = [vargin{j} '_' xname];
  eval(['MeanIRF=oo_.PosteriorIRF.Mean.' name,';']);

  if max(abs(MeanIRF)) > 1e-6 ,
    nplo=nplo+1;
    if mod(nplo,9)==1,
      figure('name',['Orthogonalised shocks to ',xname])
      nfig=nfig+1;
    end
    subplot(3,3,nplo)
    plot([1 length(MeanIRF)],[0 0],'-r','linewidth',0.5);
    %       for k = 1:9
    %         plot(1:options_.irf,DistribIRF(:,k),'-g','linewidth',0.5)
    %       end
    hold on,
    plot(1:length(MeanIRF),MeanIRF(:)./xscale,'-k','linewidth',1)
    xlim([1 length(MeanIRF)]);
    hold off
    title(vargin{j},'interpreter','none')
  end
    if (mod(nplo,9)==0 | j==length(vargin)) & nplo,
      nplo=0;
    end
end
