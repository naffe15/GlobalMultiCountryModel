% estimated folder
clear all
close all

str=[pwd '/'];
% mydrive = '/Volumes/MACFIS';
% mydrive = 'Z:';
if exist('str')~=1
if ismac
    str = '/Volumes/MACFIS/Global_Estimated_Model/GM/GM3-DE/';
else
    str ='Z:\Global_Estimated_Model\GM\GM3-DE';    
end
end
cd(str)
info=dir(str);
info1=info([info.isdir]==1);
foldername=char(info1.name);
fileID=fopen('estimated_folders_list.txt','w');
j=0;
for i=1:size(foldername,1)
    cd(deblank(foldername(i,:)))
    if isdir(['gemc' filesep 'metropolis'])
        cd (['gemc' filesep 'metropolis'])
        if exist('gemc_mh_history_0.mat')==2 && exist('metropolis.log')==2
            load gemc_mh_history_0 record
            if record.MCMCConcludedSuccessfully==1
                j=j+1;
                fprintf(fileID,'%s\n',strcat(str,foldername(i,:)));
                folder_list{j,1}=strcat(str,foldername(i,:));
                info2(j) = info1(i);
            end
        end
    end
    cd(str)
end
fclose(fileID);