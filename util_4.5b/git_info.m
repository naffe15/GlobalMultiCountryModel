function [info, git_data] = git_info(mypath, flag_notes)

[junk, git_status0]=system(['git -C ' mypath ' status -b -uno']);
if ismac
    [junk, git_log0]=system(['TERM=ansi git -C ' mypath ' log -1']);
    git_log0=regexprep(git_log0,'[[\d;]*[mK]','');
    git_log0=regexprep(git_log0,'[ ]+',' ');
else
    [junk, git_log0]=system(['git -C ' mypath ' log -1']);
end

skipline()
disp('COMMON FOLDER:')
fprintf('%s\n',mypath);
skipline()
disp('GIT STATUS:')
fprintf('%s',git_status0);
skipline()
disp('GIT LOG:')
fprintf('%s',git_log0);
skipline()
disp('END GIT INFO')
skipline()

if nargin<2
    flag_notes=0;
end

fid = fopen('git_info.txt','w');
fprintf(fid,'COMMON FOLDER:\n');
fprintf(fid,'%s\n\n', mypath);
fprintf(fid,'GIT STATUS:\n');
fprintf(fid,'%s\n',git_status0);
fprintf(fid,'GIT LOG:\n');
fprintf(fid,'%s',git_log0);
fclose(fid);

[junk, git_status]=system(['git -C ' mypath ' status -b --porcelain -uno']);
git_branch = regexp(git_status,'## .*','match','dotexceptnewline');
git_branch=regexprep(git_branch,'#','');
git_modified = regexp(git_status,'M .*','match','dotexceptnewline');

if ismac
    [junk, git_log]=system(['TERM=ansi git -C ' mypath ' log --pretty=format:"%h ::: %an, %ar : %s" -1']);
    git_log=regexprep(git_log,'[[\d;]*[mK]','');
    git_log=regexprep(git_log,'[ ]+',' ');
else
    [junk, git_log]=system(['git -C ' mypath ' log --pretty=format:"%h ::: %an, %ar : %s" -1']);
    
end
out = regexp(git_log,' ::: ','split');
git_log_commit = deblank(out{1});
git_log_commit = regexprep(git_log_commit,'[\s]','');
out = regexp(out{2},', ','split');
git_log_commit_author = out{1};
out = regexp(out{2},' : ','split');
git_log_commit_date = out{1};
git_log_commit_txt = out{2};

git_data.status=git_status;
git_data.branch=git_branch;
git_data.modified=git_modified;
git_data.log=git_log;
git_data.log_commit=git_log_commit;
git_data.log_commit_author=git_log_commit_author;
git_data.log_commit_date=git_log_commit_date;
git_data.log_commit_txt=git_log_commit_txt;
info = 0;

if flag_notes
    system(['git -C ' mypath ' fetch origin refs/notes/commits:refs/notes/origin/commits']);
    system(['git -C ' mypath ' notes merge -v origin/commits -s cat_sort_uniq']);
    %     git fetch origin refs/notes/commits:refs/notes/origin/commits
    %     git notes merge -v origin/commits -s union
    %     git notes merge -v origin/commits -s cat_sort_uniq
    if ~isempty(git_modified)
        
        skipline()
        disp('GIT warning:: your git repository is modified!')
        disp('GIT warning:: it will be impossible to trace back the commit used for this execution!')
        
        %     reply = input('Do you want to COMMIT before running the instance? Y/N [Y]:','s');
        %     if isempty(reply)
        %        reply = 'Y';
        %     end
        %     if strcmpi(reply,'Y')
        %            info = 1;
        %         error('EXECUTION STOPPED!:: COMMIT your git repo')
        %            return
        %     end
        
        %     if info
        %         skipline()
        %         disp('EXECUTION STOPPED!:: COMMIT your git repo')
        %         break
        %     else
        %         skipline()
        %         disp('EXECUTION continues:: it will be impossible to trace back the commit used for this execution')
        %         continue
        %     end
        
    else
        system(['git -C ' mypath ' notes append -m "' pwd '"']);
        system(['git -C ' mypath ' push origin refs/notes/commits']);
        %     regexp(pwd,filesep,'split')
    end
end