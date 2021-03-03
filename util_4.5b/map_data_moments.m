function loss = map_data_moments(OutputDirectoryName, Model, DynareOptions, DynareResults, EstimatedParameters, BayesInfo)

% Copyright (C) 2014 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

fname_ = Model.fname;
pnames = Model.param_names(EstimatedParameters.param_vals(:,1),:);
pvalue_ks = DynareOptions.opt_gsa.pvalue_ks;
indx_moment = [];
DynareOptions.nodisplay = 1;
init = ~DynareOptions.opt_gsa.load_stab;

options_mcf.pvalue_ks = DynareOptions.opt_gsa.pvalue_ks;
options_mcf.pvalue_corr = DynareOptions.opt_gsa.pvalue_corr;
options_mcf.alpha2 = DynareOptions.opt_gsa.alpha2_stab;
options_mcf.param_names = pnames;
options_mcf.fname_ = fname_;
options_mcf.OutputDirectoryName = OutputDirectoryName;

skipline()
disp('Sensitivity analysis for matching data moments')

if DynareOptions.opt_gsa.ppost,
    filetoload=dir([Model.dname filesep 'metropolis' filesep fname_ '_param_irf*.mat']);
    lpmat=[];
    for j=1:length(filetoload),
        load([Model.dname filesep 'metropolis' filesep fname_ '_param_irf',int2str(j),'.mat'])
        lpmat = [lpmat; stock];
        clear stock
    end
    type = 'post';
else
    if DynareOptions.opt_gsa.pprior
        filetoload=[OutputDirectoryName '/' fname_ '_prior'];
        load(filetoload,'lpmat','lpmat0','istable','iunstable','iindeterm','iwrong' ,'infox')
        lpmat = [lpmat0 lpmat];
        type = 'prior';
    else
        filetoload=[OutputDirectoryName '/' fname_ '_mc'];
        load(filetoload,'lpmat','lpmat0','istable','iunstable','iindeterm','iwrong' ,'infox')
        lpmat = [lpmat0 lpmat];
        type = 'mc';
    end
end
[Nsam, np] = size(lpmat);
npar = size(pnames,1);
nshock = np - npar;

variables_ = cellstr(DynareOptions.varobs);
nvar = length(variables_);
ar = DynareOptions.ar;
moment_data_restrictions = set_moment_data_restrictions(variables_,ar);
nbr_moment_restrictions = size(moment_data_restrictions,1);
DynareOptions.endogenous_prior_restrictions.moment = moment_data_restrictions;
DynareOptions.endogenous_prior_restrictions.irf = [];

if init
    
    mat_moment=cell(nbr_moment_restrictions,1);
    for ij=1:nbr_moment_restrictions,
        mat_moment{ij}=NaN(Nsam,length(DynareOptions.endogenous_prior_restrictions.moment{ij,3}));
    end
    
    irestrictions = [1:Nsam];
    h = dyn_waitbar(0,'Please wait...');
    for j=1:Nsam,
        Model = set_all_parameters(lpmat(j,:)',EstimatedParameters,Model);
        [Tt,Rr,SteadyState,info,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults);
        if info(1)==0,
            [info, info_irf, info_moment, data_irf, data_moment]=endogenous_prior_restrictions(Tt,Rr,Model,DynareOptions,DynareResults);
            if ~isempty(info_moment)
                for ij=1:nbr_moment_restrictions,
                    mat_moment{ij}(j,:)=data_moment{ij}(:,2)';
                end
                indx_moment(j,:)=info_moment(:,1);
            end
        else
            irestrictions(j)=0;
        end
        dyn_waitbar(j/Nsam,h,['MC iteration ',int2str(j),'/',int2str(Nsam)])
    end
    dyn_waitbar_close(h);
    
    irestrictions=irestrictions(find(irestrictions));
    xmat=lpmat(irestrictions,:);
    skipline()
    endo_prior_restrictions=DynareOptions.endogenous_prior_restrictions;
    save([OutputDirectoryName,filesep,fname_,'_',type,'_moment_match'],'xmat','mat_moment','irestrictions','indx_moment','endo_prior_restrictions');
else
    load([OutputDirectoryName,filesep,fname_,'_',type,'_moment_match'],'xmat','mat_moment','irestrictions','indx_moment','endo_prior_restrictions');
end

if ~isempty(indx_moment)
    skipline()
    disp('Deleting old MOMENT calibration plots ...')
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_match*.eps']);
    for j=1:length(a),
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_match*.fig']);
    for j=1:length(a),
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_match*.pdf']);
    for j=1:length(a),
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_match_restrictions.eps']);
    for j=1:length(a),
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_match_restrictions.fig']);
    for j=1:length(a),
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    a=dir([OutputDirectoryName,filesep,fname_,'_',type,'_moment_match_restrictions.pdf']);
    for j=1:length(a),
        delete([OutputDirectoryName,filesep,a(j).name]);
    end
    disp('done !')
    skipline()
    
    options_mcf.param_names = char(BayesInfo.name);
    all_moment_couples = cellstr([char(endo_prior_restrictions.moment(:,1)) char(endo_prior_restrictions.moment(:,2))]);
    moment_couples = unique(all_moment_couples);
    nbr_moment_couples = size(moment_couples,1);
    plot_indx = NaN(nbr_moment_couples,1);
    time_matrix=cell(nbr_moment_couples,1);
    indx_moment_matrix=zeros(length(irestrictions),nbr_moment_couples);
    moment_matrix=cell(nbr_moment_couples,1);
    moment_mean=cell(nbr_moment_couples,1);
    moment_median=cell(nbr_moment_couples,1);
    moment_var=cell(nbr_moment_couples,1);
    moment_HPD=cell(nbr_moment_couples,1);
    moment_distrib=cell(nbr_moment_couples,1);
    % For single legend search which has maximum nbr of restrictions
    maxijv=0;
    for ij=1:nbr_moment_restrictions
        if length(endo_prior_restrictions.moment{ij,3})>maxijv
            maxij=ij;maxijv=length(endo_prior_restrictions.moment{ij,3});
        end
        plot_indx(ij) = find(strcmp(moment_couples,all_moment_couples(ij,:)));
        time_matrix{plot_indx(ij)} = [time_matrix{plot_indx(ij)} endo_prior_restrictions.moment{ij,3}];
    end
    iplot_indx = ones(size(plot_indx));
    
    indx_moment = indx_moment(irestrictions,:);
    
    for ij=1:nbr_moment_restrictions,
        mat_moment{ij}=mat_moment{ij}(irestrictions,:);
        moment_matrix{plot_indx(ij)} = [moment_matrix{plot_indx(ij)} mat_moment{ij}];
        indx_moment_matrix(:,plot_indx(ij)) = indx_moment_matrix(:,plot_indx(ij)) + indx_moment(:,ij);
        for ik=1:size(mat_moment{ij},2),
            [Mean,Median,Var,HPD,Distrib] = ...
                posterior_moments(mat_moment{ij}(:,ik),0,DynareOptions.mh_conf_sig);
            moment_mean{plot_indx(ij)} = [moment_mean{plot_indx(ij)}; Mean];
            moment_median{plot_indx(ij)} = [moment_median{plot_indx(ij)}; Median];
            moment_var{plot_indx(ij)} = [moment_var{plot_indx(ij)}; Var];
            moment_HPD{plot_indx(ij)} = [moment_HPD{plot_indx(ij)}; HPD];
            moment_distrib{plot_indx(ij)} = [moment_distrib{plot_indx(ij)}; Distrib'];
        end
        leg = num2str(endo_prior_restrictions.moment{ij,3}(1));
        aleg = num2str(endo_prior_restrictions.moment{ij,3}(1));
        if size(mat_moment{ij},2)>1,
            leg = [leg,':' ,num2str(endo_prior_restrictions.moment{ij,3}(end))];
            aleg = [aleg,'_' ,num2str(endo_prior_restrictions.moment{ij,3}(end))];
            iplot_indx(ij)=0;
        end
        indx1 = find(indx_moment(:,ij)==0);
        indx2 = find(indx_moment(:,ij)~=0);
        atitle0=[endo_prior_restrictions.moment{ij,1},' vs ',endo_prior_restrictions.moment{ij,2}, '(', leg,')'];
        fprintf(['%4.1f%% of the ',type,' support matches MOMENT ',atitle0,' inside [%4.1f, %4.1f]\n'],length(indx1)/length(irestrictions)*100,endo_prior_restrictions.moment{ij,4})
        % aname=[type '_moment_match_',int2str(ij)];
        aname=[type '_moment_match_',endo_prior_restrictions.moment{ij,1},'_vs_',endo_prior_restrictions.moment{ij,2},'_',aleg];
        atitle=[type ' MOMENT Match: Parameter(s) driving ',endo_prior_restrictions.moment{ij,1},' vs ',endo_prior_restrictions.moment{ij,2}, '(', leg,')'];
        options_mcf.amcf_name = aname;
        options_mcf.amcf_title = atitle;
        options_mcf.beha_title = 'moment match';
        options_mcf.nobeha_title = 'NO moment match';
        options_mcf.title = atitle0;
        if ~isempty(indx1) && ~isempty(indx2)
            mcf_analysis(xmat, indx1, indx2, options_mcf, DynareOptions);
        end
        
        %         [proba, dproba] = stab_map_1(xmat, indx1, indx2, aname, 0);
        %         indplot=find(proba<pvalue_ks);
        %         if ~isempty(indplot)
        %             stab_map_1(xmat, indx1, indx2, aname, 1, indplot, OutputDirectoryName,[],atitle);
        %         end
    end
    
    % plot variances
    nrow=ceil(sqrt(nvar));
    ncol=nrow;
    if nrow*(nrow-1)>nvar,
        ncol=nrow-1;
    end
    hVAR=dyn_figure(DynareOptions.nodisplay,'name',[type ' evaluation of moment restrictions: VARIANCE']);
    for j=1:nvar,
        iset1 = strmatch(variables_{j},endo_prior_restrictions.moment(:,1));
        iset2 = strmatch(variables_{j},endo_prior_restrictions.moment(:,2));
        iset = intersect(iset1,iset2);
        iVAR(j) = iset(find(cell2mat(endo_prior_restrictions.moment(iset,3))==0));
        
        figure(hVAR),
        subplot(nrow,ncol,j),
        hc = cumplot(mat_moment{iVAR(j)});
        set(hc,'color','k','linewidth',2)
        hold all,
        %     hist(mat_moment{ij}),
        a=axis;
        x1val=max(endo_prior_restrictions.moment{iVAR(j),4}(1),a(1));
        x2val=min(endo_prior_restrictions.moment{iVAR(j),4}(2),a(2));
        hp = patch([x1val x2val x2val x1val],a([3 3 4 4]),'b');
        set(hp,'FaceAlpha', 0.5)
        hold off,
        title(variables_{j},'interpreter','none'),
        
        iACF = iset(find(cell2mat(endo_prior_restrictions.moment(iset,3))));
        iALL = union(iset1,iset2);
        iCCF = iALL(~ismember(iALL,iset));
        
        hCCF=dyn_figure(DynareOptions.nodisplay,'name',[type ' evaluation of ACF/CCF for variable ' variables_{j}]);
        for iv = 1:nvar,
            subplot(nrow,ncol,iv),
            if iv~=j
                isetX = iCCF(strmatch(variables_{iv},endo_prior_restrictions.moment(iCCF,1)));
                isign=1;
                if isempty(isetX)
                    isetX = iCCF(strmatch(variables_{iv},endo_prior_restrictions.moment(iCCF,2)));
                    isign=-1;
                end
                tmp_matrix=[];
                tmp=[];
                tmp_hpd = nan(ar*2+1,2);
                tmp_distrib = nan(ar*2+1,9);
                tmp_median=[];
                tmp_mean=[];
                tmp_var=[];
                for ilag=-ar:ar,
                    itmp = isetX(find(cell2mat(endo_prior_restrictions.moment(isetX,3))*isign==ilag));
                    tmp_matrix = [tmp_matrix mat_moment{itmp}];
                    [tmp_mean(ilag+ar+1),tmp_median(ilag+ar+1),tmp_var(ilag+ar+1),tmp_hpd(ilag+ar+1,:),tmp_distrib(ilag+ar+1,:)] = ...
                        posterior_moments(mat_moment{itmp},0,DynareOptions.mh_conf_sig);
                    tmp(ilag+ar+1,:) = endo_prior_restrictions.moment{itmp,4};
                end
                if ar
                    plot([-ar:ar],[max(tmp_matrix)' min(tmp_matrix)'],'k--','linewidth',2)
                    hold on,
                    plot([-ar:ar],tmp_median,'k','linewidth',2)
                    plot([-ar:ar],[tmp_distrib],'k-')
                    a=axis;
                    tmp(isinf(tmp(:,1)),1)=a(3);
                    tmp(isinf(tmp(:,2)),2)=a(4);
                    hp = patch([[-ar:ar] [ar:-1:-ar]],[tmp(:,1); tmp(end:-1:1,2)],'b');
                    set(hp,'FaceAlpha',[0.5])
                    plot(a(1:2),[0 0],'r')
                    hold off,
                    axis(a)
                    box on,
                    set(gca,'xtick',[-ar:ar])
                    title([variables_{j},'(t) ',variables_{iv} '(t-lag)'],'interpreter','none'),
                else
                    hc = cumplot(tmp_matrix);
                    set(hc,'color','k','linewidth',2)
                    hold all,
                    %     hist(mat_moment{ij}),
                    a=axis;
                    x1val=max(tmp(1),a(1));
                    x2val=min(tmp(2),a(2));
                    hp = patch([x1val x2val x2val x1val],a([3 3 4 4]),'b');
                    set(hp,'FaceAlpha', 0.5)
                    hold off,
                    title([variables_{j},' ',variables_{iv} ],'interpreter','none'),
                end
            end
            if iv==j && ar
                tmp_matrix=[];
                tmp=[1 1];
                tmp_hpd = nan(ar,2);
                tmp_distrib = nan(ar,9);
                tmp_median=[];
                tmp_mean=[];
                tmp_var=[];
                for ilag=1:ar,
                    itmp = iACF(find(cell2mat(endo_prior_restrictions.moment(iACF,3))==-ilag));
                    tmp_matrix = [tmp_matrix mat_moment{itmp}];
                    [tmp_mean(ilag),tmp_median(ilag),tmp_var(ilag),tmp_hpd(ilag,:),tmp_distrib(ilag,:)] = ...
                        posterior_moments(mat_moment{itmp},0,DynareOptions.mh_conf_sig);
                    tmp(ilag+1,:) = endo_prior_restrictions.moment{itmp,4};
                end
                plot([0:ar],[[1;max(tmp_matrix)'] [1;min(tmp_matrix)']],'k--','linewidth',2)
                hold on,
                plot([0:ar],[1 tmp_median],'k','linewidth',2)
                plot([0:ar],[ones(1,9); tmp_distrib],'k-')
                a=axis;
                tmp(isinf(tmp(:,1)),1)=a(3);
                tmp(isinf(tmp(:,2)),2)=a(4);
                hp = patch([[0:ar] [ar:-1:0]],[tmp(:,1); tmp(end:-1:1,2)],'b');
                set(hp,'FaceAlpha',[0.5])
                plot(a(1:2),[0 0],'r')
                hold off,
                axis(a)
                box on,
                set(gca,'xtick',[0:ar])
                title([variables_{j},'(t) ',variables_{j} '(t-lag)'],'interpreter','none'),
            elseif iv==j
                title([variables_{j},' ',variables_{j}],'interpreter','none'),
            end
            
        end
        dyn_saveas(hCCF,[OutputDirectoryName,filesep,fname_,'_',type,'_moment_match_restrictions_CCF_',variables_{j}],DynareOptions.nodisplay,DynareOptions.graph_format);
        
        % MCF of all ACF CCF with logical AND
        indx_set = sum(indx_moment(:,union(iACF,iCCF)),2);
        indx1 = find(indx_set==0);
        indx2 = find(indx_set~=0);
        fprintf(['%4.1f%% of the ',type,' support matches moments for ',variables_{j},'\n'],length(indx1)/length(irestrictions)*100)
        % aname=[type '_moment_match_',int2str(ij)];
        aname=[type '_moment_match_',variables_{j},'_ALL'];
        atitle=[type ' MOMENT match: Parameter(s) driving ALL moments for ',variables_{j}];
        options_mcf.amcf_name = aname;
        options_mcf.amcf_title = atitle;
        options_mcf.beha_title = 'moment match';
        options_mcf.nobeha_title = 'NO moment match';
        options_mcf.title = atitle0;
        if ~isempty(indx1) && ~isempty(indx2)
            mcf_analysis(xmat, indx1, indx2, options_mcf, DynareOptions);
        end
        
    end
    % MCF of all variances with logical AND
    indx_set = sum(indx_moment(:,iVAR),2);
    indx1 = find(indx_set==0);
    indx2 = find(indx_set~=0);
    fprintf(['%4.1f%% of the ',type,' support matches ALL variances\n'],length(indx1)/length(irestrictions)*100)
    % aname=[type '_moment_match_',int2str(ij)];
    aname=[type '_moment_match_VAR_ALL'];
    atitle=[type ' MOMENT match: Parameter(s) driving ALL variances'];
    options_mcf.amcf_name = aname;
    options_mcf.amcf_title = atitle;
    options_mcf.beha_title = 'VAR match';
    options_mcf.nobeha_title = 'NO VAR match';
    options_mcf.title = atitle0;
    if ~isempty(indx1) && ~isempty(indx2)
        mcf_analysis(xmat, indx1, indx2, options_mcf, DynareOptions);
    end
    dyn_saveas(hVAR,[OutputDirectoryName,filesep,fname_,'_',type,'_moment_match_restrictions_VAR'],DynareOptions.nodisplay,DynareOptions.graph_format);
    
    if ar,
        for ij=1:nbr_moment_couples,
            itmp = min(find(plot_indx==ij));

            % MCF of the couples with logical AND
            indx1 = find(indx_moment_matrix(:,ij)==0);
            indx2 = find(indx_moment_matrix(:,ij)~=0);
            leg = num2str(time_matrix{ij}(1));
            leg = [leg '...' num2str(time_matrix{ij}(end))];
            aleg = 'ALL';
            atitle0=[endo_prior_restrictions.moment{itmp,1},' vs ',endo_prior_restrictions.moment{itmp,2}, '(', leg,')'];
            fprintf(['%4.1f%% of the ',type,' support matches MOMENT restrictions ',atitle0,'\n'],length(indx1)/length(irestrictions)*100)
            % aname=[type '_moment_calib_',int2str(ij)];
            aname=[type '_moment_match_',endo_prior_restrictions.moment{itmp,1},'_vs_',endo_prior_restrictions.moment{itmp,2},'_',aleg];
            atitle=[type ' MOMENT match: Parameter(s) driving ',endo_prior_restrictions.moment{itmp,1},' vs ',endo_prior_restrictions.moment{itmp,2}, '(', leg,')'];
            options_mcf.amcf_name = aname;
            options_mcf.amcf_title = atitle;
            options_mcf.beha_title = 'moment match';
            options_mcf.nobeha_title = 'NO moment match';
            options_mcf.title = atitle0;
            if ~isempty(indx1) && ~isempty(indx2)
                mcf_analysis(xmat, indx1, indx2, options_mcf, DynareOptions);
            end
        end
    end
    
    skipline()
end
return

