function clean_dynare_proj(varargin)
% clean_dynare_proj(varargin)

if nargin==0,
    disp('clean_dynare_proj(NAME_of_DYNARE_Project)')
    return
end

for j=1:length(varargin),
    disp(['Deleting ',varargin{j},' folder ...'])
    [SUCCESS,MESSAGE,MESSAGEID] = rmdir(varargin{j},'s');
    disp('Done!')
    disp(['Deleting ',varargin{j},'*.mat files ...'])
    delete([varargin{j},'*.mat'])
    disp('Done!')
    disp(['Deleting ',varargin{j},'*.eps files ...'])
    delete([varargin{j},'*.eps'])
    disp('Done!')
    disp(['Deleting ',varargin{j},'*.fig files ...'])
    delete([varargin{j},'*.fig'])
    disp('Done!')
    disp(['Deleting ',varargin{j},'*.pdf files ...'])
    delete([varargin{j},'*.pdf'])
    disp('Done!')
    disp(['Deleting ',varargin{j},'.m _dynamic.m _static.m files ...'])
    delete([varargin{j},'.m'])
    delete([varargin{j},'_dynamic.m'])
    delete([varargin{j},'_static.m'])
    disp('Done!')
end