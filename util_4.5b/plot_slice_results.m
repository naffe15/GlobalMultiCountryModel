ifig=0;
for indx=1:nt,
    if mod(indx,9)==1,
        ifig=ifig+1;
        figure(ifig),
        iplo=0;
    end
    iplo=iplo+1;
    [post_mean, post_median, post_var, hpd_interval, post_deciles, ...
        density] = posterior_moments(XSIM(:,indx),1,options_.mh_conf_sig);
    [x,f,abscissa,dens,binf,bsup] = draw_prior_density(indx,bayestopt_);
    subplot(3,3,iplo),
    plot(abscissa,dens,'LineWidth',2,'Color',[0.7 0.7 0.7])
    hold all, plot(density(:,1),density(:,2), 'k','LineWidth',2)
    yl=get(gca,'ylim');
    set(gca,'ylim', [0 yl(2)]),
    set(gca,'xlim', [bayestopt_.lb(indx) bayestopt_.ub(indx)]),
    plot(xparam1([indx indx]),[0 yl(2)],'g')
    title(bayestopt_.name{indx},'interpreter','none')
end
