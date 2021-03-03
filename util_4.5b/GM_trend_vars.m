function GM_trend_vars(countries, T0, datafile)

load(datafile)

T=T0;

eval(['divP=exp(GP0*((1-tbase):(length(T)-tbase)))'';']);
eval(['divPOP=exp(GPOP0*(1:length(T))'');']);
for c=1:length(countries)
    co=countries{c};
    eval(['divR=meanYOBS_',co,'.*exp(GYTREND0*(1:length(T))'');']);
    eval(['divR_',co,'= divR;']);
    eval(['divPOP_',co,'= divPOP*meanPOP_',co,';']);
    if isempty(strmatch('RoW',countries{c}))
        eval(['divN_',co,'= divPOP.*meanPOP_',co,'.*meanHPERE_',co,'.*meanACTR_',co,'.*meanPARTR_',co,'.*meanL_',co,';']);
        eval(['divPC=exp(GPC0_',co,'*((1-tbase):(length(T)-tbase)))'';']);
        eval(['divPC_',co,'=divPC;']);
        eval(['divPG=exp(GPG0_',co,'*((1-tbase):(length(T)-tbase)))'';']);
        eval(['divPG_',co,'=divPG;']);
        eval(['divPI=exp(GPITOT0_',co,'*((1-tbase):(length(T)-tbase)))'';']);
        eval(['divPI_',co,'=divPI;']);
        eval(['divTFP_',co,'= (divR./divN_',co,').^(0.65)./(divP./divPI).^(1-0.65);']);
    end
end

save trendfile div*