dname=ls('*.mod');
dname=dname(1:end-4);
delete('*.eps')
delete('*.pdf')
delete('*.asv')
delete('*.fig')
delete('*.tex')
delete('*.log*')
flag=rmdir(dname,'s');

