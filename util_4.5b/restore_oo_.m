fnam = fieldnames(oo_sd);
for k=1:length(fnam)
    oo_.(fnam{k}) = oo_sd.(fnam{k});
end