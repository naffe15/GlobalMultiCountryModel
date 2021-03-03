function dynare2git(origDIR, gitDIR, commit)
%function dynare2git(origDIR, gitDIR, commit)
%
% copy of dynare project files to GIT repository
% commit = 1 optional automatic commit, message is just name of directory

if nargin<3 
    commit=0;
end

fname = ls(gitDIR);
for j=1:size(fname,1),
    if ~isdir([gitDIR,filesep,fname(j,:)]) && ~isempty(dir([origDIR,filesep,fname(j,:)])),
        copyfile([origDIR,filesep,fname(j,:)], gitDIR)
    end
end
if commit
    currDIR=pwd;
    indx=strfind(origDIR,filesep);
    msg = origDIR(indx(end)+1:end);
    cd(gitDIR);
    status = system(['"C:\Program Files (x86)"\Git\bin\git commit -a -m ', msg]);
    if status
        system(['"C:\Program Files"\Git\bin\git commit -a -m ', msg]);
    end
    cd(currDIR);
end
