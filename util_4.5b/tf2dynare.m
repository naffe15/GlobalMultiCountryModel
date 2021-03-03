function tf2dynare(DENDEN, NUMNUM, fname, nfwd, nbkwd)

[nout, nlag, nin] = size(NUMNUM);
if nargin<4,
  nfwd = floor(length(DENDEN)/2);
  nbkwd = nfwd;
end

fid = fopen([fname,'_emul.mod'],'w');
str = 'var ';
for j=1:nout, str=[str,'y',int2str(j),' ']; end
str=[str, '\n'];
fprintf(fid,str);
str = '';
for j=1:nin, str=[str,'u',int2str(j),' ']; end
str=[str, ';\n\n'];
fprintf(fid,str);
str = 'varexo ';
for j=1:nin, str=[str,'e',int2str(j),' ']; end
str=[str, ';\n\n'];
fprintf(fid,str);
str = 'parameters ';
for j=1:nbkwd, str=[str,'a',int2str(j),' ']; end
for j=1:nfwd, str=[str,'c',int2str(j),' ']; end
str=[str, '\n'];
fprintf(fid,str);
for i=1:nout,
  str = '';
  for k=1:nin,
    for j=1:nlag, str=[str,'b_',int2str(i),'_',int2str(k),'_',int2str(j-1),' ']; end
    str=[str, '\n'];
  end
  fprintf(fid,str);
end
str = '';
for k=1:nin,
  str=[str, 'rho',int2str(k),' '];
end
  fprintf(fid,str);
fprintf(fid,';\n\n\n');

for j=1:nbkwd,
  fprintf(fid,['a',int2str(j),'=%27.18e;\n'],DENDEN(j+nfwd+1));
end  
for j=1:nfwd,
  fprintf(fid,['c',int2str(j),'=%27.18e;\n'],DENDEN(nfwd+1-j));
end  
for j=1:nin,
  fprintf(fid,['rho',int2str(j),'=0;\n']);
end  

for i=1:nout,
  for k=1:nin,
    for j=1:nlag,     
      fprintf(fid,['b_',int2str(i),'_',int2str(k),'_',int2str(j-1),' =%27.18e;\n'],NUMNUM(i,j,k));
    end
  end
end
fprintf(fid,'\n\n');
fprintf(fid,'model(linear);\n');
for j=1:nout,
  str = ['y',int2str(j)];
  for i=1:nfwd,
    str = [str, ' + c',int2str(i),'*y',int2str(j),'(+',int2str(i),')'];
  end
  fprintf(fid,str);
  str = ['\n    '];
  for i=1:nbkwd,
    str = [str, ' + a',int2str(i),'*y',int2str(j),'(-',int2str(i),')'];
  end
  fprintf(fid,[str,' =\n']);
  for k=1:nin,
    str = ['    '];
    str = [str, ' + b_',int2str(j),'_',int2str(k),'_0*u',int2str(k)];
    for i=1:nlag-1,     
      str = [str, ' + b_',int2str(j),'_',int2str(k),'_',int2str(i),'*u',int2str(k),'(-',int2str(i),')'];
    end
    fprintf(fid,[str,'\n']);
  end
  fprintf(fid,';\n');
end
for k=1:nin,
  fprintf(fid,['u',int2str(k),' = rho',int2str(k),'*u',int2str(k),'(-1) + e',int2str(k),';\n']);
end
fprintf(fid,'end;\n\n');
fprintf(fid,'steady;\n');
fprintf(fid,'check;\n');
fprintf(fid,'M_.Sigma_e=eye(M_.exo_nbr);\n');
fprintf(fid,'stoch_simul(periods=2000, irf=200) y1;\n');


fclose(fid);
