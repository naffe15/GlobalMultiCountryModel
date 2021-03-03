function [paths_avg,res,bounds,res1]=probam(paths,prc,neigv)
%
% Input: paths(i,:) is the i^th path, 
% 'prc' is a row vector containing percents (20=20%),
% 'neigv' is the third dimension of res(:,:,:), see below.
%
% Output: paths_avg is the average of paths (row vector). 
% For each i, paths(i,:) can be expressed as 
% paths_avg+sum_j{gamma(i,j)*eigvec(:,j)'}, where
% eigvec(:,j)' is the j^th eigenvector.
% For a given percent prc(l), let maki_lk denote the prc(l) percentile
% of gamma(:,k). Then res(l,:,k)=paths_avg+maki_lm*eigvec(:,m)', where m is
% the index of the eigenvector corresponding to the k^th largest eigenvalue. 
% bounds is simple (see below). 
% If m is the index of the eigenvector corresponding to the largest
% eigenvalue, and gamma(i,m) is the element of gamma(:,m) corresponding
% to the prc(l) percentile of gamma(:,m), then res1(l,:)=paths(i,:).
% Other details are in probam.pdf 

[npaths,lpaths]=size(paths);
nprc=length(prc);     

bounds=zeros(2,lpaths);
prcb=[prc(1) prc(nprc)];
for i=1:lpaths
[prctemp1,prctemp2]=prctil(paths(:,i),prcb);
bounds(:,i)=prctemp1';
end

w=cov(paths);
[eigvec,eigval]=eig(w);
eigvals=diag(eigval);
paths_avg=zeros(1,lpaths);
for i=1:npaths
    paths_avg=paths_avg+paths(i,:);
end
paths_avg=paths_avg/npaths;
gamma=zeros(npaths,lpaths);
for i=1:npaths
    gamma(i,:)=(eigvec'*(paths(i,:)-paths_avg)')';
end
[eigsrt,eigord]=sort(eigvals);

res=zeros(nprc,lpaths,neigv);
for i=1:neigv
  [prctemp1,prctemp2]=prctil(gamma(:,eigord(lpaths-i+1)),prc);
  res(:,:,i)=(eigvec(:,eigord(lpaths-i+1))*prctemp1)';
  for j=1:nprc
  res(j,:,i)=res(j,:,i)+paths_avg;
  end
end

[prctemp1,prctemp2]=prctil(gamma(:,eigord(lpaths)),prc);
res1=paths(prctemp2(:),:);

function [prcx,prcxk]=prctil(prcy,prcz)
%prcx returns prcz percentiles of vector prcy, and prcxk
%is such that prcx=prcy(prcxk)
ly=length(prcy);
lz=length(prcz);
prcx=prcz;
prcxk=prcz;
[temp,ix]=sort(prcy); %temp=prcy(ix)
prcz=(ly/100)*prcz;
for i=1:lz
    if prcz(i)<ly/2
        prcz(i)=ceil(prcz(i));
    else
        prcz(i)=ly-ceil(ly-prcz(i))+1;
    end
    prcx(i)=temp(prcz(i)); %prcx=prcy(ix(prcz))
    prcxk(i)=ix(prcz(i));
end

