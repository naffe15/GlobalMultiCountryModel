function vardec_groups(Vi, Vinf, obsname, ex_names_, leg, file_name)
global M_

if nargin<6 || isempty(file_name)
    file_name='';
end

nobs = size(Vi,1);
nfcast = size(Vi,3);
ngroups = size(ex_names_,1);
Vi_group = zeros(nobs,ngroups+1,nfcast);
Vinf_group  = zeros(nobs,ngroups+1);
Vtmp=Vi;
Vtmp2=Vinf;
for j=1:ngroups,
    indx = find(ismember(cellstr(M_.exo_names),ex_names_(j,:)));
    Vi_group(:,j,:)=sum(Vi(:,indx,:),2);
    Vinf_group(:,j)=sum(Vinf(:,indx),2);
    Vtmp(:,indx,:)=0;
    Vtmp2(:,indx)=0;
end
Vi_group(:,ngroups+1,:)=sum(Vtmp,2);
Vinf_group(:,ngroups+1)=sum(Vtmp2,2);

file_name  = [M_.fname,'_vardec_',int2str(ngroups),'shocks',file_name,'.xls'];
delete(file_name);
iprint=[1, 4, 8, 16, 40, 100, 200, 400];
xlstxt=[ leg(1:end-1,1) ex_names_];
if ~ismac
    xlswrite(file_name,xlstxt,'legend');
else
    xlswrite_MACOS(file_name,xlstxt,'legend');
end
    for j=1:length(iprint),
    
    ip=iprint(j);
    if nfcast>=ip,      
        
        xlstxt=[ [cellstr(' ') leg(:,1)'];[obsname num2cell(Vi_group(:,:,ip))]];
        if ~ismac
            xlswrite(file_name,xlstxt,[int2str(ip),'-step']);
        else
            xlswrite_MACOS(file_name,xlstxt,[int2str(ip),'-step']);
        end
    end
end

title=['VARIANCE DECOMPOSITION (inf-step ahead, in percent)'];
        
xlstxt=[[cellstr(' ') leg(:,1)']; [obsname num2cell(Vinf_group)]];
if ~ismac
    xlswrite(file_name,xlstxt,['Inf-step']);
else
    xlswrite_MACOS(file_name,xlstxt,['Inf-step']);
end

