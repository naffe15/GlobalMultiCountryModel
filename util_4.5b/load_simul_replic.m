function yy = load_simul_replic(var_list_, fname)
global oo_ options_ M_

[nr nc]=size(oo_.endo_simul);
nc = nc-M_.maximum_lag;

if nargin<2 | isempty(fname),
  fname=[M_.fname,'_simul'];
end


for i=1:size(var_list_,1),
fh=fopen(fname);
  iv=strmatch(deblank(var_list_(i,:)), M_.endo_names, 'exact');
for  j=1:options_.replic,
  A=fread(fh,[nr nc],'float64');
  yy(j,:)=A(iv,:);
  
end
fclose(fh)
end