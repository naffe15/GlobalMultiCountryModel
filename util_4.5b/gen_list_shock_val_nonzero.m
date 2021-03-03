function gen_list_shock_val_nonzero(M_,fname,OutDir),

% generates command  lines that assigns values currently stored in M_.Sigma_e (e.g.
% estimated values) to model shocks


if nargin<2,
    disp('shocks;')
for j=1:M_.exo_nbr, 
    if M_.Sigma_e(j,j)
        disp(['var ',deblank(M_.exo_names(j,:)),';']), 
        disp(['stderr ',num2str(sqrt(M_.Sigma_e(j,j)),'%21.14g'),';']), 
    end
end
    disp('end;')

else
if nargin<3,
    OutDir='.';
else
    OutDir = CheckPath(OutDir,M_.fname);
end
    fid = fopen([OutDir filesep M_.fname,'_',fname],'w');
    fprintf(fid,'shocks;\n');
for j=1:M_.exo_nbr, 
    if M_.Sigma_e(j,j)
        fprintf(fid,['var ',deblank(M_.exo_names(j,:)),';\n']);
        fprintf(fid,'stderr %21.14g;\n',sqrt(M_.Sigma_e(j,j))); 
    end
end
    fprintf(fid,'end;\n');
fclose(fid);
end    