function [bigVAR, yy1, yy2, x00]=vecm2(xx, xvec, beta, nex, nk, nlag)
global opts_vecm_

a=ones(1,size(xx,2));
b=a*0;
b(nex)=ones(length(nex),1);
c=a-b;
nendo=find(c);
nvar=size(xx,2);
nvec=size(xvec,2);

% for j=1:nvar;
%     x0(j) =  mean(xx(:,j));
%     xx(:,j)=xx(:,j)-x0(j);
% end
% 
% if nvec,
%     x0vec=mean(xvec);
%     xvec=xvec-ones(length(xvec),1)*mean(xvec);
% else
%     x0vec=[];    
% end
% 
% exogenous sub-block
ind_ex_=1;
xxx=ones(size(xx,1)-nlag,1);
for j=1:nlag,
    xxx = [xxx xx(j:end-nlag+j-1,nex)];
    ind_ex_ = [ind_ex_ 1+nex+nvar*(j-1)];
end
yyy = xx(nlag+1:end,nex) ;
bigVAR=zeros(nvar,nvar*nlag+nvec+1);
for j=1:length(nex),
    %bigrho(j,:) = ((xxx'*xxx)\sum(xxx.*kron(yyy(:,j),ones(1,size(xxx,2))))')';
    bigVAR(nex(j),ind_ex_) = ((xxx'*xxx)\sum(xxx.*kron(yyy(:,j),ones(1,size(xxx,2))))')';
end
if nvec,
[j0, i0]=find(beta(:,nex));
if ~isempty(j0),   
    if length(i0)==1,
        xxxx = [xxx xvec(1:end-nlag,j0)];
        bigVAR(nex(i0),[ind_ex_ nvar*nlag+1+j0])= ((xxxx'*xxxx)\sum(xxxx.*kron(yyy(:,i0),ones(1,size(xxxx,2))))')';        
    elseif ~find(diff(i0))
        xxxx = [xxx xvec(1:end-nlag,j0)];
        bigVAR(nex(i0),[ind_ex_ nvar*nlag+1+j0])= ((xxxx'*xxxx)\sum(xxxx.*kron(yyy(:,i0),ones(1,size(xxxx,2))))')';        
    else
        ii=unique(i0);        
        for j=1:length(ii),
            jj=find(i0==ii(j));
            %xxx = [xx(1:end-1,nex) xvec(1:end-1,j0(jj))];
            xxxx = [xxx xvec(1:end-nlag,j0(jj))];
            %bigVAR(i0,[nex nvar+j0(jj)])=((xxx'*xxx)\sum(xxx.*kron(yyy(:,j),ones(1,size(xxx,2))))')';
            bigVAR(nex(ii(j)),[ind_ex_ nvar*nlag+1+j0(jj)])= ((xxxx'*xxxx)\sum(xxxx.*kron(yyy(:,ii(j)),ones(1,size(xxxx,2))))')';        
        end
    end
end
end
%bigrho = ((xxx'*xxx)\xxx'*yyy)';

% full VECM
xxx=ones(size(xx,1)-nlag,1);
for j=1:nlag,
    xxx=[xxx xx(j:end-nlag+j-1,:)];
end
if nvec,
    xxx=[xxx xvec(1:end-nlag,:)];
end
ind_vec_ = [nvar*nlag+2:nvar*nlag+nvec+1];
yyy=xx(nlag+1:end,nendo);

for j=1:length(nendo),
    bigVAR(nendo(j),:) = ((xxx'*xxx)\sum(xxx.*kron(yyy(:,j),ones(1,size(xxx,2))))')';
end
if nvec,
dumVECM=beta*bigVAR(:,1);
for j=1:nlag,
    dumVECM = [dumVECM beta*bigVAR(:,2+nvar*(j-1):nvar*j+1)];
end
dumVECM=[dumVECM  eye(size(beta,1))+beta*bigVAR(:,ind_vec_) ];

if nlag>1,
bigVAR = [bigVAR; ...
        [zeros(nvar*(nlag-1), 1) eye(nvar*(nlag-1)) zeros(nvar*(nlag-1), nvec+nvar)]; ...   % more than 1 lag rows
        dumVECM];   % equilibrium correcton raws
else
    bigVAR = [bigVAR; ...
        dumVECM];   % equilibrium correcton raws
end
end

yy=(bigVAR*xxx')';  % 1-step ahead prediction
yy=yy(:,1:(nvar+nvec));
if nvec,
    yyy=[xx(nlag+1:end,:) xvec(nlag+1:end,:)]; % data
else
    yyy=[xx(nlag+1:end,:)]; % data
end
    
   
dum0 = (eye(nlag*nvar+nvec)-bigVAR(:,2:end))\bigVAR(:,1);
x00=dum0([[1:nvar] ind_vec_-1]);
xxx=[];
for j=1:nlag,
    %xxx=[xxx [zeros(j-1,nvar); ...
    %    xx(j:end,:)] ];
    xxx=[xxx [ones(j-1,1)*x00(1:nvar)'; ...
        xx(j:end,:)] ];
end
xxx=[xxx xvec];
yy1=ones(size(xxx,1)+nk,nvar+nvec,nk);
for j=1:length(x00),
    yy1(:,j,:)=yy1(:,j,:).*x00(j);
end

for j=1:nk,
   dum=xxx;
   for i=1:j,
       dum=[ones(size(xx,1),1) dum];
       dum=(bigVAR*dum')';
   end
   yy1(j+1:size(xxx,1)+j,:, j)=dum(:,[[1:nvar] ind_vec_-1]);
end

for j=1:opts_vecm_.forecast,
    dum=xxx(end,:);
   for i=1:j, 
       dum=[1 dum];
       dum=(bigVAR*dum(end,:)')';
   end
   yy2(j,:)=dum([[1:nvar] ind_vec_-1]);
end

if max(sort(abs(eig(bigVAR(:,2:end)))))>1,
    disp('WARNING: the unconstrained mean VECM is unstable!!!')
end