function []=WriteParameters2Excel(M_, oo_, bayestopt_, xlsfilename)

if nargin<4,
    xlsfilename=[M_.fname,'_Parameters.xls'];
end

if ~exist('oo_') || ~exist('M_') 
    error(['quest3hlmr_results.mat'])
end


EstimatedParamsNames=fieldnames(oo_.posterior_mode.parameters);
EstimatedShocksNames=fieldnames(oo_.posterior_mode.shocks_std);
ParamsTexNames=M_.param_names_tex;
ShocksTexNames=M_.exo_names_tex;
ParamsLongNames=M_.param_names_long;
ShocksLongNames=M_.exo_names_long;


NbrOfEstimatedParams=size(EstimatedParamsNames,1);
NbrOfEstimatedShocks=size(EstimatedShocksNames,1);

AllParamsNames=M_.param_names;

CheckIfEstimated=ones(size(AllParamsNames,1),1);

for jj=1:size(AllParamsNames,1)
    TempParName=deblank({AllParamsNames(jj,:)});
    for kk=1:size(EstimatedParamsNames,1)
        TempEstimParName=EstimatedParamsNames(kk,:);
        
        if strcmp(TempParName,TempEstimParName)
            CheckIfEstimated(jj,1)=0;
        end
    end
end



AllShocksNames=M_.exo_names;

CheckIfShockEstimated=ones(size(AllShocksNames,1),1);

for jj=1:size(AllShocksNames,1)
    TempParName={deblank(AllShocksNames(jj,:))};
    for kk=1:NbrOfEstimatedShocks
        TempEstimParName=EstimatedShocksNames(kk,:);
        
        if strcmp(TempParName,TempEstimParName)
            CheckIfShockEstimated(jj,1)=0;
        end
    end
end


CalibratedShocksIndx=find(CheckIfShockEstimated==1);
EstimatedShocksIndx=find(CheckIfShockEstimated==0);

CalibratedParamsIndx=find(CheckIfEstimated==1);
EstimatedParamsIndx=find(CheckIfEstimated==0);

CalibratedParams={};
CalibratedParamsValues=M_.params(CalibratedParamsIndx);
CalibratedParamsNames=M_.param_names(CalibratedParamsIndx,:);
CalibratedParamsTexNames=M_.param_names_tex(CalibratedParamsIndx,:);
CalibratedParamsLongNames=M_.param_names_long(CalibratedParamsIndx,:);
%EstimatedParamsTexNames=M_.param_names_tex(EstimatedParamsIndx,:);
EstimatedParamsTexNames=[];
for i=1:length(EstimatedParamsNames)
    ii=strmatch(EstimatedParamsNames{i},AllParamsNames,'exact');
    EstimatedParamsTexNames{i}=ParamsTexNames(ii,:);
    EstimatedParamsLongNames{i}=ParamsLongNames(ii,:);
end
%EstimatedShocksTexNames= M_.exo_names_tex(EstimatedShocksIndx,:);
EstimatedShocksTexNames=[];
for i=1:length(EstimatedShocksNames)
    ii=strmatch(EstimatedShocksNames{i},AllShocksNames,'exact');
    EstimatedShocksTexNames{i}=ShocksTexNames(ii,:);
    EstimatedShocksLongNames{i}=ShocksLongNames(ii,:);
end

CalibratedParams(1,1)={'Calibrated Parameters'};
CalibratedParams(1,2)={'Latex name'};
CalibratedParams(1,3)={'calibrated value'};
CalibratedParams(1,4)={'Long name'};
for jj=1:size(CalibratedParamsIndx,1)
    CalibratedParams(jj+1,3)=num2cell(CalibratedParamsValues(jj,1));
    CalibratedParams(jj+1,1)={deblank(CalibratedParamsNames(jj,:))};
    %     CalibratedParams(jj+1,2)={CalibratedParamsTexNames(jj,:)};
    TempTexName=deblank(CalibratedParamsTexNames(jj,:));
    TempTexName=TempTexName(2:end-1);
    CalibratedParams(jj+1,2)={TempTexName};
    CalibratedParams(jj+1,4)={CalibratedParamsLongNames(jj,:)};
end




if nargin<3 || isempty('bayestopt_')
    load([M_.fname '\prior\definition.mat'])
end

PriorType=bayestopt_.pshape;
PriorMean=bayestopt_.p1;
PriorStd=bayestopt_.p2;

ShocksPriorType=PriorType(1:NbrOfEstimatedShocks);
ShocksPriorMean=PriorMean(1:NbrOfEstimatedShocks);
ShocksPriorStd=PriorStd(1:NbrOfEstimatedShocks);

ParamsPriorType=PriorType(NbrOfEstimatedShocks+1:end);
ParamsPriorMean=PriorMean(NbrOfEstimatedShocks+1:end);
ParamsPriorStd=PriorStd(NbrOfEstimatedShocks+1:end);


EstimatedParams={};
EstimatedShocks={};


EstimatedParams(1,1)={'Estimated Parameters'};
EstimatedParams(1,2)={'Latex name'};
EstimatedParams(1,3)={'Prior type'};
EstimatedParams(1,4)={'Prior Mean'};
EstimatedParams(1,5)={'Prior Std'};
EstimatedParams(1,6)={'Posterior Mode'};
EstimatedParams(1,7)={'Posterior Std'};
EstimatedParams(1,8)={'Long name'};

for jj=1:NbrOfEstimatedParams
    EstimatedParams(jj+1,1)=EstimatedParamsNames(jj,:);
    TempTexName=deblank(EstimatedParamsTexNames{jj});
    TempTexName=TempTexName(2:end-1);
    
    EstimatedParams(jj+1,2)={TempTexName};
    EstimatedParams(jj+1,8)={EstimatedParamsLongNames{jj}};
    
    if ParamsPriorType(jj)==1;
        EstimatedParams(jj+1,3)={'Beta'};
    elseif ParamsPriorType(jj)==2
        EstimatedParams(jj+1,3)={'Gamma'};
    end
    EstimatedParams(jj+1,4)=num2cell(ParamsPriorMean(jj,1));
    EstimatedParams(jj+1,5)=num2cell(ParamsPriorStd(jj,1));
    EstimatedParams(jj+1,6)=num2cell(getfield(oo_.posterior_mode.parameters, EstimatedParamsNames{jj,:}));
    EstimatedParams(jj+1,7)=num2cell(getfield(oo_.posterior_std_at_mode.parameters, EstimatedParamsNames{jj,:}));
end

NbrOfShocks=M_.exo_nbr;
CalibratedShocksIndx = setxor(1:NbrOfShocks, EstimatedShocksIndx);
NbrOfCalibratedShocks = length(CalibratedShocksIndx);
CalibratedShocks={};
CalibratedShocks(1,1)={'Calibrated Shocks'};
CalibratedShocks(1,2)={'Latex name'};
CalibratedShocks(1,3)={'Value'};
CalibratedShocks(1,4)={'Long name'};
for jj=1:NbrOfCalibratedShocks
    CalibratedShocks(jj+1,1)={M_.exo_names(CalibratedShocksIndx(jj),:)};
    
    TempTexName=deblank(M_.exo_names_tex(CalibratedShocksIndx(jj),:));
    TempTexName=TempTexName(2:end-1);    
    CalibratedShocks(jj+1,2)={TempTexName};
    
    TempLongName=M_.exo_names_long(CalibratedShocksIndx(jj),:);
    CalibratedShocks(jj+1,4)={TempLongName};
    
    TempValue=sqrt(M_.Sigma_e(CalibratedShocksIndx(jj),CalibratedShocksIndx(jj)));
    CalibratedShocks(jj+1,3)=num2cell(TempValue);
end

EstimatedShocks(1,1)={'Estimated Shocks'};
EstimatedShocks(1,2)={'Latex name'};
EstimatedShocks(1,3)={'Prior type'};
EstimatedShocks(1,4)={'Prior Mean'};
EstimatedShocks(1,5)={'Prior Std'};
EstimatedShocks(1,6)={'Posterior Mode'};
EstimatedShocks(1,7)={'Posterior Std'};
EstimatedShocks(1,8)={'Long name'};

for jj=1:NbrOfEstimatedShocks
    EstimatedShocks(jj+1,1)=EstimatedShocksNames(jj,:);
    
    TempTexName=deblank(EstimatedShocksTexNames{jj});
    TempTexName=TempTexName(2:end-1);    
    EstimatedShocks(jj+1,2)={TempTexName};

    TempLongName=(EstimatedShocksLongNames{jj});    
    EstimatedShocks(jj+1,8)={TempLongName};
    
    if ShocksPriorType(jj)==1;
        EstimatedShocks(jj+1,3)={'Beta'};
    elseif ShocksPriorType(jj)==2
        EstimatedShocks(jj+1,3)={'Gamma'};
    end
    EstimatedShocks(jj+1,4)=num2cell(ShocksPriorMean(jj,1));
    EstimatedShocks(jj+1,5)=num2cell(ShocksPriorStd(jj,1));
    EstimatedShocks(jj+1,6)=num2cell(getfield(oo_.posterior_mode.shocks_std, EstimatedShocksNames{jj,:}));
    EstimatedShocks(jj+1,7)=num2cell(getfield(oo_.posterior_std_at_mode.shocks_std, EstimatedShocksNames{jj,:}));
end


if ~ismac
    xlswrite(xlsfilename,CalibratedParams,'Calibrated Params')
    xlswrite(xlsfilename,CalibratedShocks,'Calibrated Shocks')
    xlswrite(xlsfilename,EstimatedParams,'Estimated Params')
    xlswrite(xlsfilename,EstimatedShocks,'Estimated Shocks')
else
%     xlswrite_MACOS(xlsfilename,CalibratedParams,'Calibrated Params')
%     xlswrite_MACOS(xlsfilename,CalibratedShocks,'Calibrated Shocks')
%     xlswrite_MACOS(xlsfilename,EstimatedParams,'Estimated Params')
%     xlswrite_MACOS(xlsfilename,EstimatedShocks,'Estimated Shocks')
    writetable(array2table(CalibratedParams), [xlsfilename ], 'Sheet', 'Calibrated Params'); 
    writetable(array2table(CalibratedShocks), [xlsfilename ], 'Sheet', 'Calibrated Shocks'); 
    writetable(array2table(EstimatedParams), [xlsfilename ], 'Sheet', 'Estimated Params'); 
    writetable(array2table(EstimatedShocks), [xlsfilename ], 'Sheet', 'Estimated Shocks'); 
 
end