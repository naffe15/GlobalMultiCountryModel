function equation_tags_latex(M_, filename)

if nargin==1 || isempty(filename),
    filename = [M_.fname,'_dynamic.tex'];
end
fid  = fopen(filename);
aa=fscanf(fid,'%c',Inf);
fclose(fid);

ggg=strfind(aa,'\begin{dmath}');
bb=aa(1:ggg(1)-1);
candidates = M_.equations_tags(strmatch('name',M_.equations_tags(:,2),'exact'),:);
indices = cell2mat( candidates(:,1) );
for j=1:length(ggg),
    indx = find(indices==j);
    if ~isempty(indx)
        bb = [bb candidates{indx,3} sprintf('\n')];
    end
    if j<length(ggg),
       bb = [bb aa(ggg(j):(ggg(j+1)-1))];
    else
       bb = [bb aa(ggg(j):end)];
    end
end
% aa = strrep(aa, '{dmath}', '{dmath}');
% aa = strrep(aa, '\usepackage{fullpage}', '\usepackage{fullpage,breqn}');
f1=fopen(filename,'w');
% f1=fopen('quest3hlmr_formulae.tex','w');
fprintf(f1,'%c',bb);
fclose(f1);