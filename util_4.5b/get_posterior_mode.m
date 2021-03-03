function x=get_posterior_mode
global oo_

fn=fieldnames(oo_.posterior_density.shocks_std);
for j=1:size(fn,1),
  xx=getfield(oo_.posterior_density.shocks_std,fn{j});
  [dum, i]=max(xx(:,2));
  x(j,1)=xx(i,1);
end
offset=length(x);
fn=fieldnames(oo_.posterior_density);
for j=1:size(fn,1)-1,
  xx=getfield(oo_.posterior_density,fn{j});
  [dum, i]=max(xx(:,2));
  x(j+offset,1)=xx(i,1);
end

