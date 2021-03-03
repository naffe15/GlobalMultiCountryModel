function assign_options(opts)
   
option_names = fieldnames(opts);
for j=1:length(option_names)
    assignin('caller',option_names{j}, opts.(option_names{j}));
end

end
