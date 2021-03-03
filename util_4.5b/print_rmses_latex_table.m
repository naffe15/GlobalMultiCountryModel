function print_rmses_latex_table(obs0,rmse,r2);

global options_ bayestopt_ M_ oo_


nobs=length(obs0);

texname = M_.endo_names_tex;
lgy_ = M_.endo_names;

mfys=[];
for j=1:length(obs0),
    dum = strmatch(obs0{j},lgy_,'exact');
    nametex(j,:)=texname(dum,:);
    mfys = [mfys dum];
    if j==1,
        lgobs_ = obs0{j};
    else
        lgobs_ = str2mat(lgobs_,obs0{j});
    end
end





fid=fopen('RmsesTable.TeX','w');


fprintf(fid,'%s\n','{\tiny');

fprintf(fid,'%s\n','\begin{table}');

fprintf(fid,'%s\n','\centering');

fprintf(fid,'%s\n','\begin{tabular}{l|lcc}');

fprintf(fid,'%s\n','\hline\hline \\ ');


fprintf(fid,'%s\n','  & RMSE& $r^2$. \\');
fprintf(fid,'%s\n','\hline\\ ');

for i=1:nobs
    fprintf(fid,'%s',  '$' );
    fprintf(fid,'%s',  nametex(i,:) );
    fprintf(fid,'%s',  '$' );
    fprintf(fid,'%s','&');
    fprintf(fid,'%s', num2str(rmse(i)) );
    fprintf(fid,'%s','&');
    fprintf(fid,'%s', num2str(r2(i)) );
    fprintf(fid,'%s\n','\\ ');
    
end

fprintf(fid,'%s\n','\hline\hline');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n',['\caption{RMSE and $r^2$}']);
fprintf(fid,'%s\n','\end{table}');
fprintf(fid,'%s\n','}');







fclose(fid);

% open('RmsesTable.TeX')
