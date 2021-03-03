function wrap_eps_figs_for_latex(varargin)

delete('wrap_eps_figs.TeX')


fidTeX = fopen(['wrap_eps_figs.TeX'],'w+');
if ~ismac
    list_of_eps = ls('*.eps');
else
    list_of_eps = dir('*.eps');
end
if ~ismac
for j=1:size(list_of_eps,1),
    tempfilename = list_of_eps(j,:);
    [a,b,c]=fileparts(tempfilename);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,['\\includegraphics[width=0.80\\textwidth] {' b '} \n']);
    fprintf(fidTeX,['\\caption{\\protect\\url{' deblank(tempfilename) '}}\n']);
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,' \n');
end

for i=1:nargin
    if ~isempty(dir([varargin{i} '/*.eps']))
        list_of_eps = ls([varargin{i} '/*.eps']);
        for j=1:size(list_of_eps,1),
            tempfilename = list_of_eps(j,:);
            [a,b,c]=fileparts(tempfilename);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,['\\includegraphics[width=0.80\\textwidth] {' varargin{i} '/' b '} \n']);
            fprintf(fidTeX,['\\caption{\\protect\\url{' varargin{i} '/' deblank(tempfilename) '}}\n']);
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
        end
    end
end
else
    for j=1:size(list_of_eps,1),
    tempfilename = list_of_eps(j).name;
    [a,b,c]=fileparts(tempfilename);
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,['\\includegraphics[width=0.80\\textwidth] {' b '} \n']);
    fprintf(fidTeX,['\\caption{\\protect\\url{' deblank(tempfilename) '}}\n']);
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,' \n');
end

for i=1:nargin
    if ~isempty(dir([varargin{i} '/*.eps']))
        list_of_eps = dir([varargin{i} '/*.eps']);
        for j=1:size(list_of_eps,1),
            tempfilename = list_of_eps(j).name;
            [a,b,c]=fileparts(tempfilename);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,['\\includegraphics[width=0.80\\textwidth] {' varargin{i} '/' b '} \n']);
            fprintf(fidTeX,['\\caption{\\protect\\url{' varargin{i} '/' deblank(tempfilename) '}}\n']);
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
        end
    end
end
end

fclose(fidTeX);



