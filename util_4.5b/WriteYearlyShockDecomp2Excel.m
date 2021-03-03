function out=WriteYearlyShockDecomp2Excel(M_,MatFileList)

if  nargin <2
    MatFileList=[M_.fname,'_YoY_ShockDecomp.mat'];
end
MatFileList = dir(MatFileList);

out=[];
for jf=1:length(MatFileList),
    MatFileName = MatFileList(jf).name;
tempstruct=load(MatFileName);

ExistingFields=fieldnames(tempstruct);




for jj=1:size(ExistingFields,1)
    a2=getfield(tempstruct,ExistingFields{jj,:});
        XlsFileName = [a2.TagName '.xls'];
        if exist(XlsFileName,'file') && jj==1,
            delete(XlsFileName);
        end
    Decomp=a2.Decomp';
    Init_Decomp=a2.Init_Decomp';
    Smooth=a2.Smooth';
    Init_Smooth=a2.Init_Smooth';
    
    d0={};
    
    NbrOfShocks=size(Decomp,2);
    LengthOfVar=size(Decomp,1);
    Shocks=a2.Shocks;
    if size(Shocks,2)>1 || ischar(Shocks{1}),
        old_interface=1;
    elseif iscell(Shocks{1})
        old_interface=0;
        ngrps=size(Shocks,1);
        ncc = 0;
        for j=1:ngrps,
            ncc = max(ncc,length(Shocks{j}));
        end
        Shocks_=cell([ngrps ncc]);
        for j=1:ngrps,
            Shocks_(j,1:length(Shocks{j})) = Shocks{j};
        end
        Shocks=Shocks_;
        clear Shocks_;
    end
    
    sizeshocks=size(Shocks);
    Names=a2.Names(:,1)';
    InitialTime=a2.InitialTime;
    EndTime=InitialTime+fix(LengthOfVar/4)+(LengthOfVar-4*fix(LengthOfVar/4))/4;
    TimeSpan=InitialTime:0.25:EndTime;
    TimeSpan=TimeSpan(1:end-1)';
    
    d0(1,1)={'Init Decomposition'};
    
    
    d0(1,NbrOfShocks+2)={'Init Smoot Var'};
    d0(1,NbrOfShocks+3)={'No Init Decomposition'};
    d0(1,2:NbrOfShocks+1)=Names;
    d0(1,NbrOfShocks+4:2*NbrOfShocks+3)=Names;
    d0(2:LengthOfVar+1,2:NbrOfShocks+1)=num2cell(Init_Decomp);
    d0(2:LengthOfVar+1,NbrOfShocks+4:2*NbrOfShocks+3)=num2cell(Decomp);
    
    d0(2:LengthOfVar+1,NbrOfShocks+2)=num2cell(Init_Smooth);
    d0(2:LengthOfVar+1,2*NbrOfShocks+4)=num2cell(Smooth);
    LastRow=size(d0,1);
    d0(LastRow+2,1)={'Legend:'};
    d0(LastRow+2,2)={'Shocks include:'};
    d0(LastRow+3:LastRow+3+sizeshocks(1),1)=Names';
    d0(LastRow+3:LastRow+2+sizeshocks(1),2:sizeshocks(2)+1)=Shocks;
    d0(1,end)={'Init Smooth Var'};
    d0(2:length(TimeSpan)+1,1)=num2cell(TimeSpan);
    d0(2,NbrOfShocks+3)=num2cell(a2.InitialTime);
        warning off
        if ~ismac
            [STATUS,MESSAGE] = xlswrite(XlsFileName,d0, a2.SheetName);
        else
            [STATUS,MESSAGE] = xlswrite_MACOS(XlsFileName,d0, a2.SheetName);
        end
        warning on
    clear d0
end
end