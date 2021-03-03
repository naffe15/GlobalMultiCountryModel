function [bigVAR, yy1, yy2, xm]=vecm0(xx, xvec, beta, nex, nk)
global opts_vecm_

a=ones(1,size(xx,2));
b=a*0;
b(nex)=ones(length(nex),1);
c=a-b;
nendo=find(c);
nvar=size(xx,2);
nvec=size(xvec,2);

for j=1:nvar;
   if opts_vecm_.thmean
       x0(j) =  opts_vecm_.y0(j); 
   else
       x0(j) =  mean(xx(:,j)); 
   end
   xx(:,j)=xx(:,j)-x0(j);
end

if nvec,
    x0vec=mean(xvec);
    xvec=xvec-ones(length(xvec),1)*mean(xvec);
else
    x0vec=[];    
end

bigVAR=zeros(nvar,nvar+nvec);
% exogenous sub-block
xxx = [xx(1:end-1,nex)];
yyy = xx(2:end,nex) ;
bigVAR(nex,nex) = ((xxx'*xxx)\xxx'*yyy)';
if nvec,
[j0, i0]=find(beta(:,nex));
if ~isempty(j0),   
    if length(i0)==1,
        xxx = [xx(1:end-1,nex) xvec(1:end-1,j0)];
        bigVAR(nex(i0),[nex nvar+j0])= ((xxx'*xxx)\sum(xxx.*kron(yyy(:,(i0)),ones(1,size(xxx,2))))')';        
    elseif ~find(diff(i0))
        xxx = [xx(1:end-1,nex) xvec(1:end-1,j0)];
        bigVAR(nex(i0),[nex nvar+j0])= ((xxx'*xxx)\sum(xxx.*kron(yyy(:,(i0)),ones(1,size(xxx,2))))')';        
    else
        ii=unique(i0);        
        for j=1:length(ii),
            jj=find(i0==ii(j));
            xxx = [xx(1:end-1,nex) xvec(1:end-1,j0(jj))];
            bigVAR(nex(ii(j)),[nex nvar+j0(jj)])=((xxx'*xxx)\sum(xxx.*kron(yyy(:,(ii(j))),ones(1,size(xxx,2))))')';
        end
    end
end
end

% full VECM
if nvec,
    xxx=[xx(1:end-1,:) xvec(1:end-1,:)];
else
    xxx=[xx(1:end-1,:)];
end
yyy=xx(2:end,nendo);


for j=1:size(yyy,2),
    bigVAR(nendo(j),:) = ((xxx'*xxx)\sum(xxx.*kron(yyy(:,j),ones(1,size(xxx,2))))')';
end
if nvec,
bigVAR = [bigVAR; ...
    beta*bigVAR(:,1:nvar),  (eye(size(beta,1))+beta*bigVAR(:,nvar+1:end) )];   % equilibrium correcton raws
end

%xxx=[xx(1:end-1,:) xvec(1:end-1,:)];
yy=(bigVAR*xxx')';  % 1-step ahead prediction
if nvec,
    yyy=[xx(2:end,:) xvec(2:end,:)]; % data
else
    yyy=xx(2:end,:); % data
end
   
   
xxx=[xx xvec];
x00=[x0 x0vec];
yy1=ones(size(xxx,1)+nk,size(xxx,2),nk);
for j=1:length(x00),
    yy1(:,j,:)=yy1(:,j,:).*x00(j);
end
for j=1:nk,
   yy1(j+1:size(xxx,1)+j,:, j)=(bigVAR^j*xxx')' + ones(size(xxx,1),1)*[x0 x0vec];
end

for j=1:1:opts_vecm_.forecast,
   yy2(j,:)=(bigVAR^j*xxx(end,:)')'+[x0 x0vec];
end

xm=[x0 x0vec];

if max(sort(abs(eig(bigVAR))))>1,
    disp('WARNING: the VECM is unstable!!!')
end
