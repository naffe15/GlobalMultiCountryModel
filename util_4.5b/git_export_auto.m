for j=1:length(info2)
    if ismac
        [~,git_hash]=system(['TERM=ansi git -C ../git log --grep=' info2(j).name ':: --format="%H"']);
    else
        [~,git_hash]=system(['git -C ../git log --grep=' info2(j).name ':: --format="%H"']);
    end
    if length(git_hash)>8
        warning(['git_export_auto:: commit for ' info2(j).name ':: already exists!'])
        disp('git_export_auto:: skip commit.')
        continue
    end
    
    copyfile([folder_list{j,1} filesep 'gemc.dyn'],'.')
    [status,result] = system('git add gemc.dyn');
    if status
        error(result)
    end
    [status,result] = system(['git commit -m "' info2(j).name '::" --allow-empty']);
    if status
        error(result)
    end
    
end