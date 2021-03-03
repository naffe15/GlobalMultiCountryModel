function gen_list_param_val(M_,fname,OutDir),

% generates command  lines that assigns values currently stored in M_ (e.g.
% estimated values) to model params

if nargin<2,
for j=1:M_.param_nbr, 
    disp([M_.param_names(j,:),' = ',num2str(M_.params(j),'%21.14g'),';']), 
end

else
if nargin<3,
    OutDir='.';
else
    OutDir = CheckPath(OutDir,M_.fname);
end
    fid = fopen([OutDir filesep M_.fname,'_',fname],'w');
for j=1:M_.param_nbr, 
    fprintf(fid,[M_.param_names(j,:),' = %21.14g;\n'],M_.params(j)); 
end
fclose(fid);
end    