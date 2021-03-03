function M_ = setLatexName(M_,vname,texname);

ncendo = size(        M_.endo_names_tex ,2);
ncexo = size(        M_.exo_names_tex ,2);

for j=1:length(vname),
    jendo = strmatch(vname{j},M_.endo_names);
    if ~isempty(jendo),
        M_.endo_names_tex(jendo,1:ncendo) = ' ';
        M_.endo_names_tex(jendo,1:length(texname{j})) = texname{j};
    end
    jexo = strmatch(vname{j},M_.exo_names);
    if ~isempty(jexo),
        M_.exo_names_tex(jexo,1:ncexo) = ' ';
        M_.exo_names_tex(jexo,1:length(texname{j})) = texname{j};
    end
end



