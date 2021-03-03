papertype=1;
if papertype
dirname='ImbalancesPaper';
else
 dirname='ACESPaper';   
end
if isempty(dir(dirname))
    mkdir(dirname)
end

smoothfrcslist_0=ls('*smooth_and*.fig');

for jj=1:size(smoothfrcslist_0,1)
    copyfile(smoothfrcslist_0(jj,:), [dirname,'\'])
end


smoothfrcslist_assa01=ls('quest3hlmr/assa01/*smooth_and*.fig');
for jj=1:size(smoothfrcslist_assa01,1)
    oldname=strcat('quest3hlmr\assa01\',smoothfrcslist_assa01(jj,:));
    
    newname=strcat(dirname,'\HigerRiskPremium_',smoothfrcslist_assa01(jj,:));
    copyfile(oldname, newname)
end


smoothfrcslist_assaTL=ls('quest3hlmr/assaTL/*smooth_and*.fig');
for jj=1:size(smoothfrcslist_assaTL,1)
   oldname=strcat('quest3hlmr\assaTL\',smoothfrcslist_assa01(jj,:));
    
    newname=strcat(dirname,'\FiscalConsolidation_',smoothfrcslist_assaTL(jj,:));
    copyfile(oldname, newname)
end

smoothfrcslist_assaTFP=ls('quest3hlmr/assaTFP/*smooth_and*.fig');
for jj=1:size(smoothfrcslist_assaTFP,1)
   oldname=strcat('quest3hlmr\assaTFP\',smoothfrcslist_assa01(jj,:));
    
    newname=strcat(dirname,'\LowerGrowth_',smoothfrcslist_assaTFP(jj,:));
    copyfile(oldname, newname)
end


listshocks=ls('*Unobserved*Varia*.fig');

for jj=1:size(listshocks,1)
    
    copyfile(listshocks(jj,:), [dirname,'\'])
end

listshockdecs=ls('*9Shocks*Init.fig');

for jj=1:size(listshockdecs,1)
    
    copyfile(listshockdecs(jj,:), [dirname,'\'])
end


copyfile('Decomposition_New.fig', [dirname,'\'])
