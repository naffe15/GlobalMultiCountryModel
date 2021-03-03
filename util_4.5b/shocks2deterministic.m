function shocks2deterministic(etahat, M_, options_, file_name)

if nargin<3
    file_name = 'det_shocks.dyn';
end
fid = fopen(file_name,'w+') ;
 
fprintf(fid,'shocks;\n');

for jx = 1:M_.exo_nbr
    fprintf(fid,'var %s;\n', M_.exo_names(jx,:));
    
    fprintf(fid,'periods ');
    fprintf(fid,' %2i',1:options_.forecast);
    fprintf(fid,';\n');
    
    fprintf(fid,'values ');
    fprintf(fid,' (%12.8e)',etahat(jx,end-options_.forecast+1:end));
    fprintf(fid,';\n');
end
fprintf(fid,'end;\n');

fclose(fid);
