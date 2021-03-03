function print_cond_var_latex_table(matrix,ex_names_,names_shock,rowLabels)

global options_

filename = 'var_dec.tex';
textsize = 'footnotesize';
width   = size(ex_names_,1)+1;
height  = size(matrix, 1);
matrix  = round(matrix,2);
nvars   = size(rowLabels,2);
steps   = options_.conditional_variance_decomposition(1:end-1)-1;



fid = fopen(filename, 'w');
fprintf(fid, '\\begin{table}\r\n');
fprintf(fid, '\\setlength{\\tabcolsep}{1.25em}\r\n');
if(~isempty(textsize))
    fprintf(fid, '\\begin{%s}', textsize);
end
fprintf(fid, '\\begin{tabularx}{\\linewidth}\r\n');
colLabels = names_shock;

% Row label
fprintf(fid, '{l');

for i=1:width
   % fprintf(fid, '%c|', alignment);
    fprintf(fid, 'X');        
end
fprintf(fid, '}\r\n');

fprintf(fid, '\\toprule\r\n');


% Col Label
if(~isempty(colLabels))
    if(~isempty(rowLabels))
        fprintf(fid, '&');
    end
    for w=1:width-1
        %fprintf(fid, '\\textbf{%s}&', colLabels{w});
        fprintf(fid, '%s &', colLabels{w});
    end
        fprintf(fid, '%s\\\\\\midrule\r\n', colLabels{width});
end
var_count=0;
% Rows 
for h=1:height
 
    if h == 1
        fprintf(fid, '\\multicolumn{%1.0f}{l}{\\textit{Conditional %1.0f-step ahead}}\\\\\r\n', width(1)+1,steps(1));
    end
    
    if h == nvars+1
        fprintf(fid, '\\multicolumn{%1.0f}{X}{}\\\\\r\n', width(1)+1);
        fprintf(fid, '\\multicolumn{%1.0f}{l}{\\textit{Conditional %1.0f-step ahead}}\\\\\r\n', width(1)+1,steps(2));
    end
    
    if h == nvars*2+1
        fprintf(fid, '\\multicolumn{%1.0f}{X}{}\\\\\r\n', width(1));
        fprintf(fid, '\\multicolumn{%1.0f}{l}{\\textit{Unconditional variance}}\\\\\r\n', width(1)+1);
    end
    
    if(~isempty(rowLabels))
        if var_count == size(rowLabels,2)
            var_count = 0;
        end
        var_count = var_count +1;
        fprintf(fid, '%s&', rowLabels{var_count});
    
    end
    for w=1:width-1        
        field = sprintf('%1.1f',matrix(h, w));
        fprintf(fid, '%s&', field);
    %fprintf(fid, '%s&', matrix{h, w});
    end

field = sprintf('%1.1f',matrix(h,width(1)));
fprintf(fid, '%s\\\\\r\n',field);



end
fprintf(fid, '\\bottomrule\r\n');
fprintf(fid, '\\end{tabularx}\r\n');

if(~isempty(textsize))
    fprintf(fid, '\\end{%s}', textsize);
end

fprintf(fid, '\\end{table}\r\n');
fclose(fid);
end
