function git_note_gemc()

if ismac
    [~,userid]=system('whoami');
    if isequal(deblank(userid),'iazzusr')
        excl = {'daemon','Guest', 'nobody', 'iazzusr','root'};
        [~,user_list]=system('dscl . list /Users | grep -v ''^_''');
        matches = regexp(user_list,'.*','match','dotexceptnewline');
        matches = matches(~ismember(matches,excl));
        if length(matches)==1
            userid = matches{1};
        else
            disp('there are more possible users on this mac!')
            userid = input('tell me who you are:','s');
        end
    else
        userid = 0;
    end
end


if ~isdir('../git')
    warning('git_note_gemc:: git repo does not exist!!!')
    disp('git_note_gemc:: skip note.')
    return
end
mydir=pwd;
str = regexp(mydir, filesep, 'split');
str=str{end};
if ismac
    [~,git_hash]=system(['TERM=ansi git -C ../git log --grep=' str ':: --format="%H"']);
else
    [~,git_hash]=system(['git -C ../git log --grep=' str ':: --format="%H"']);
end

if length(git_hash)<=8
    error(['git_note_gemc:: commit for ' str ':: does not exist. I cannot add a note! Run git_commit_gemc, instead.'])
    return
end

cd ../git

fname = [CreateTimeString '.txt'] ;
fid = fopen(fname,'w');
fprintf(fid,'\r\n');
fprintf(fid,'%s\r\n\r\n','# type here the note to append');
fclose(fid);
copyfile( fname, ['bkp_' fname])

if ismac
    if userid
        str1 = regexp(mydir, filesep, 'split');
        system(['sudo -u ' userid ' open -W -a TextEdit  /Volumes/MACFIS/' str1{end-3} filesep str1{end-2} filesep str1{end-1} '/git/' fname]);
    else
        [s,r]=system(['open -W -a TextEdit  ' fname]);
    end
else
    [s,r]=system(['notepad ' fname]);
end
A=dir(fname);
B=dir(['bkp_' fname]);

B.name=A.name;

if isequal(A,B)
    delete(fname)
    delete(['bkp_' fname])
    cd(mydir)
    %     error('YOU MUST EDIT THE README of THE GM SESSION!')
    return
else
    
    filetext = fileread(fname);
    
    expr = '[#].*\n';
    matches = regexp(filetext,expr,'match','dotexceptnewline');
    
    S = filetext;
    
    for k=1:length(matches)
        S = strrep(S,matches{k},'');
    end
    
    f1=fopen(fname,'w');
    fprintf(f1,'%c',S);
    fclose(f1);
    
end

[status,result] = system(['git notes append --allow-empty -F ' fname]);

if status ==0   
    f1=fopen([mydir filesep 'GM_README.txt'],'a');
    fprintf(f1,'%c',S);
    fclose(f1);
else
    warning('git_note_gemc:: note could not be added')
    warning(['git_note_gemc:: ' result])
end
delete(fname)
delete(['bkp_' fname])
cd(mydir)
