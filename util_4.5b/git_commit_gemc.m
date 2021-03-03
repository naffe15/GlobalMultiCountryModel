function git_commit_gemc(git_data)

if nargin==0
    git_data=0;
end

if ischar(git_data)
    [info, git_data] = git_info(git_data,0);
end

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
    warning('git_commit_gemc:: git repo does not exist!!!')
    disp('git_commit_gemc:: skip commit.')
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
if length(git_hash)>8
    warning(['git_commit_gemc:: commit for ' str ':: already exists!'])
    disp(['git_commit_gemc:: skip commit.'])
    return
end

copyfile gemc.dyn ../git/
cd ../git

fname = [CreateTimeString '.txt'] ;
fid = fopen(fname,'w');
fprintf(fid,'\r\n%s\r\n','# type here the README of this GM model instance');
fclose(fid);
copyfile( fname, ['bkp_' fname])

if ismac
    %     system(['gedit ' fname]);
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
    error('YOU MUST EDIT THE README of THE GM SESSION!')
    
else
    
    filetext = fileread(fname);
    
    expr = '[#].*\n';
    matches = regexp(filetext,expr,'match','dotexceptnewline');
    
    S = filetext;
    
    for k=1:length(matches)
        S = strrep(S,matches{k},'');
    end
        
end

fid = fopen(fname,'w');
fprintf(fid,'%s:: \r\n\r\n',str);
if ~isequal(git_data,0)
    fprintf(fid,'mypath info::\r\n');
    fprintf(fid,'branch:: %s\r\n',git_data.branch{:});
    if ~isempty(git_data.modified)
        fprintf(fid,'modified:: %s\r\n',git_data.modified{:});
    else
        fprintf(fid,'clean (unmodified)::\r\n');
    end
else
    fprintf(fid,'local session (no -Dmypath)::\r\n');
end
fprintf(fid,'\r\n');
fprintf(fid,'%c',S);
fclose(fid);


if ismac
    [~,git_hash]=system('TERM=ansi git log -1 --format="%H"');
else
    [~,git_hash]=system('git log -1 --format="%H"');
end

[status,result] = system('git add gemc.dyn');
if status
    [s,r]=system(['git reset --hard ' git_hash]);
    delete(fname)
    delete(['bkp_' fname])
    cd(mydir)
    error(result)
end

[status,result] = system(['git commit --allow-empty -F ' fname]);
if status
    [s,r]=system(['git reset --hard ' git_hash]);
    delete(fname)
    delete(['bkp_' fname])
    cd(mydir)
    error(result)
end

copyfile( fname, [mydir filesep 'GM_README.txt'])
delete(fname)
delete(['bkp_' fname])
cd(mydir)
