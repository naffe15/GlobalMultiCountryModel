function plot_update(varargin)
global options_ M_ oo_

ys_ = oo_.steady_state*0;
lgy_ = M_.endo_names;
fname_ = M_.fname;
texname=M_.endo_names_tex;

ifig0=0;
if strmatch('-',varargin{end})
    optn = varargin{end};
    varargin=varargin(1:end-1);
    if strmatch('-new',optn)
        ifig0=length(dir([fname_,'_UpdatedUnobserved*.fig']));
    end
else
    afig=dir([fname_,'_UpdatedUnobserved*.*']);
    for j=1:length(afig),
        delete(afig(j).name);
    end
    delete([M_.fname '_Updated_Unobserved_Plots.TeX'])
end
if exist(options_.datafile)
    instr = options_.datafile;
else
    instr = ['load ' options_.datafile];
end
try
    eval(instr);
    if ~exist('T','var'),
        temp = eval(deblank(options_.varobs(1,:)));
        T=[1:length(temp)]';
        clear temp;
    end
catch
    T=[1:options_.nobs+options_.first_obs-1];
end

fobs = options_.first_obs;
nv=length(varargin);
ttrend = [1:length(T(fobs:fobs+options_.nobs-1))]';
for ifig = 1:ceil(nv/9),
    h = dyn_figure(options_.nodisplay,'Name','Updated unobserved variables');
    for j=1+9*(ifig-1):min(9*ifig,nv),
        exponential_trend=0;
        trendvar{j} = '';
        if ~isempty(strfind(varargin{j},'('))
            i0 = strfind(varargin{j},'(');
            i1 = strfind(varargin{j},')')-1;
            trend(j) = evalin('base',varargin{j}(i0+1:i1));
            varargin{j}=varargin{j}(1:i0-1);
        elseif ~isempty(strfind(varargin{j},'['))
            i0 = strfind(varargin{j},'[');
            i1 = strfind(varargin{j},',')-1;
            if isempty(i1),
                i1 = strfind(varargin{j},']')-1;
            else
                i2 = strfind(varargin{j},']')-1;
                trendvar{j} = varargin{j}(i1+2:i2);
            end                
            trend(j) = evalin('base',varargin{j}(i0+1:i1));
            varargin{j}=varargin{j}(1:i0-1);
            exponential_trend=1;
        else
            trend(j) = 0;
        end
        subplot(3,3,j-9*(ifig-1)),
        i = strmatch(varargin{j},lgy_,'exact');
        posi(j)=i;
        if ~isempty(texname)
            if strmatch('\log',texname(i,:));
                texname(i,1:end-1)=texname(i,2:end);
            end
        end
        if exist(varargin{j},'var')
            yobs = eval([varargin{j},'(fobs:fobs+options_.nobs-1)']);
        end
        if exponential_trend==1,
            if ~isempty(trendvar{j}),
                iplus = strfind(trendvar{j},'+');
                iminus = strfind(trendvar{j},'-');
                if ~isempty(iplus),
                    techno=get_smooth(trendvar{j}(1:iplus(1)-1)) + get_smooth(trendvar{j}(iplus(1)+1:end));
                elseif ~isempty(iminus),
                    techno=get_smooth(trendvar{j}(1:iminus(1)-1)) + get_smooth(trendvar{j}(iminus(1)+1:end));
                else
                    techno = get_smooth(trendvar{j});
                end
            else
                techno=0;
            end
            ymodel = eval(['[(oo_.UpdatedVariables.',varargin{j},'(1:end)+ys_(i)).*exp(trend(j)*ttrend+techno)]']);
            if exist(varargin{j},'var') && trend(j)
                ymodel = ymodel.*mean_nan(yobs./ymodel);
            end
        else
            ymodel = eval(['[oo_.UpdatedVariables.',varargin{j},'(1:end)+ys_(i)+trend(j)*ttrend]']);
            if exist(varargin{j},'var') && trend(j)
                ymodel = ymodel+mean_nan(yobs-ymodel);
            end
        end
        plot(T(fobs:fobs+options_.nobs-1),ymodel,'-k')
        if exist(varargin{j},'var')
            hold on,
%             eval(['plot(T(fobs:fobs+options_.nobs-1),',varargin{j},'(fobs:fobs+options_.nobs-1),''-r.'')'])
            plot(T(fobs:fobs+options_.nobs-1),yobs,':r.')
        end
        
        if options_.TeX,
            title(varargin{j},'interpreter','none')
        else
            title(deblank(texname(i,:)),'interpreter','tex')
        end
        set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
        hold on, plot([T(fobs)-1 T(fobs+options_.nobs-1)+1],[get_mean(varargin{j}) get_mean(varargin{j})],':b')
    end
    dyn_saveas(gcf,[fname_,'_UpdatedUnobserved',int2str(ifig+ifig0)],options_.nodisplay,options_.graph_format)
    if options_.TeX,
        h1=1;
        numberofvars=nv;
        texname0=texname(posi,:);
        
        fidTeX = fopen([M_.fname '_Updated_Unobserved_Plots.TeX'],'a+');
        tempfilename=[fname_,'_UpdatedUnobserved',int2str(ifig+ifig0)];
        [a,b,c]=fileparts(tempfilename);
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for kk=1:min(numberofvars,9)
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(varargin{kk}),['$' deblank(texname0(h1,:)) '$']);
            h1=h1+1;
        end
        numberofvars=numberofvars-9;
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,['\\includegraphics[width=0.80\\textwidth] {' b '} \n']);
        fprintf(fidTeX,'\\caption{Updated Unobserved %s}',int2str(ifig+ifig0));
        fprintf(fidTeX,'\\label{Fig: Updated Unobserved:%s}\n',int2str(ifig+ifig0));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
        fclose(fidTeX); 
    end
end


