function info = git_notes_fetch(mypath)
% EXAMPLE
% git_notes_fetch ../../common_bea

[info1, a]=system(['git -C ' mypath ' fetch origin refs/notes/commits:refs/notes/origin/commits']);
[info2, b]=system(['git -C ' mypath ' notes merge -v origin/commits -s cat_sort_uniq']);

if nargout
    info = any([info1 info2]);
end
skipline()
disp('GIT NOTES FETCH:')
fprintf('%s',a);
fprintf('%s',b);
skipline()
