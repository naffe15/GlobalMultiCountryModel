function etahat = get_shocks_from_smoother(SmoothedShocks,M_)

fnam = fieldnames(SmoothedShocks);
for j=1:M_.exo_nbr
    etahat(j,:)=SmoothedShocks.(fnam{j});
end


