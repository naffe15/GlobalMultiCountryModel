function info = git_notes_push(mypath)
% EXAMPLE
% git_notes_push ../../common_bea

[info1, a]=system(['git -C ' mypath ' push origin refs/notes/commits']);

if nargout
    info = info1;
end
skipline()
disp('GIT NOTES PUSH:')
fprintf('%s',a);
skipline()

