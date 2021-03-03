function moment_data_restrictions = set_moment_data_restrictions(variables_, ar)


data_restrictions=data_covariances(variables_{:});
bandwidth = 0.5;
nvar = length(variables_);
icount= 0;
moment_data_restrictions = cell(nvar*(ar+1)+nvar*(nvar-1)/2*(2*ar+1),4);
for j=1:nvar,
    icount=icount+1;
    moment_data_restrictions{icount,1}=variables_{j};
    moment_data_restrictions{icount,2}=variables_{j};
    moment_data_restrictions{icount,3}=0;
    moment_data_restrictions{icount,4}=data_restrictions.var(j,j)*[1-bandwidth 1/(1-bandwidth)];
    if ar
    tmp1=getfield(data_restrictions.acf,variables_{j});
    end
    for ij=1:ar,
        icount=icount+1;
        moment_data_restrictions{icount,1}=variables_{j};
        moment_data_restrictions{icount,2}=variables_{j};
        moment_data_restrictions{icount,3}=-ij;
        moment_data_restrictions{icount,4}=data_restrictions.autocorr(j,ij)+[-tmp1(ij,3) tmp1(ij,3)];
    end
    for iv=j+1:nvar,
        tmp_data = getfield(data_restrictions.ccf,[variables_{j} '_' variables_{iv}]);
        for ij=ar:-1:-ar,
            icount=icount+1;
            moment_data_restrictions{icount,1}=variables_{j};
            moment_data_restrictions{icount,2}=variables_{iv};
            moment_data_restrictions{icount,3}=-ij;
            %moment_data_restrictions{icount,4}=sort(tmp_data(ij+ar+1,2)*[1-bandwidth 1+bandwidth]);
            moment_data_restrictions{icount,4}=tmp_data(ij+ar+1,2)+[-tmp_data(ij+ar+1,3) tmp_data(ij+ar+1,3)];
        end
    end
end
