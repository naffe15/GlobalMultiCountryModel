endonamestex=M_.endo_names_tex;


list02=[];

list01=[];
j=0;
[A,B]=dynare_resolve(M_,options_,oo_);
for i=1:size(endonamestex,1);
%   i
    if strfind(endonamestex(i,:),'AUX'),
        j=j+1;
        s1=strfind(endonamestex(i,:),'\');
        varpos1=(s1(end-1))+2;
        varpos2=(s1(end))-1;
        lagstep=str2num(endonamestex(i,varpos2+3));
%         varpos=str2num(endonamestex(i,varpos1:varpos2))+1;
%         s2(j)=varpos;
        ix=i;
        if lagstep>1,
            endotemp=endonamestex(i,:);
            endotemp(varpos2+3)='1';
            ix = strmatch(endotemp,endonamestex);
        end
        dump1=find(oo_.dr.order_var==ix);
%         dump2=find(abs(A(dump1,:))>1.e-6);
        [d1,dump2]=max(abs(A(dump1,:)));
        varpos=oo_.dr.order_var(dump2,:);
        
        endpos0=strfind(endonamestex(i,:), ' ');
        if ~isempty(endpos0)
            list01=strvcat(list01,endonamestex(i,1:endpos0-1));
        else
            list01=strvcat(list01,endonamestex(i,:));
        end
        endpos=strfind(endonamestex(varpos,:), ' ');
        
        list02=strvcat(list02, [endonamestex(varpos,1:endpos-1) '(' num2str(-lagstep) ')']);
    end,
end

list1=list01;

list2=list02;


fid  = fopen('quest3hlmr_dynamic.tex ');
aa=fscanf(fid,'%c',Inf);
fclose(fid);

for i=1:size(list1,1)
    aa = strrep(aa, deblank(list1(i,:)),deblank(list2(i,:)));
end

aa = strrep(aa, '{equation}', '{dmath}');
aa = strrep(aa, '\usepackage{fullpage}', '\usepackage{fullpage,breqn}');

f1=fopen('quest3hlmr_dynamic.tex','w');
% f1=fopen('quest3hlmr_formulae.tex','w');
fprintf(f1,'%c',aa);
fclose(f1);
