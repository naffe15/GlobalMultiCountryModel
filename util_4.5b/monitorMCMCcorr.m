function monitorMCMCcorr(j1,j2,xparam1)
global bayestopt_

if ~isnumeric(j1),
param1 = j1;
param2 = j2;
% param1 = 'SIGCE';
% param2 = 'SNLC';
j1 = strmatch(param1,bayestopt_.name,'exact');
j2 = strmatch(param2,bayestopt_.name,'exact');
end
% j1=29;
% j2=42;

fname = 'quest3hlmr';
firstchain=1;
nchain = 4;
aa=dir([fname,'/metropolis/',fname,'_mh*_blck*.mat']);
firstifile = 1;
nfiles = floor(length(aa)/nchain);
ll=[];
for j=firstchain:firstchain+nchain-1, 
    for i=firstifile:nfiles, 
    load([fname,'/metropolis/',fname,'_mh',int2str(i),'_blck',int2str(j)]); ll=[ll; logpo2]; 
end, end
nr=size(ll,1)/nchain;


xx=[]; 
for j=firstchain:firstchain+nchain-1,  
    for i=firstifile:nfiles, load([fname,'/metropolis/',fname,'_mh',int2str(i),'_blck',int2str(j)]); 
        xx=[xx; x2(:,[j1 j2])]; 
    end, 
end,

figure,
pscatter(xx,char(bayestopt_.name([j1 j2])));
figure,
ccontour(xx(:,1),xx(:,2),ll), 
xlabel(bayestopt_.name(j1))
ylabel(bayestopt_.name(j2))
if nargin==3,
hold on, plot(xparam1(j1),xparam1(j2),'r*')
end

