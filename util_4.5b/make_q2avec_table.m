function q2avec_table = make_q2avec_table(q2avec,vlist)
 q2avec_table = struct();
 xlist = {'frcst_name','plot'};
 myFields = fieldnames(q2avec(1));
 myFields=myFields(~ismember(myFields,xlist));
 NbFields = numel(myFields);

 for k=1:size(vlist,1)
     indx = strmatch(vlist{k,1},{q2avec.name});
     if isempty(indx)
         indx = strmatch(vlist{k,1},{q2avec.gname});
     end
     
     if isempty(indx)
         error(['make_q2avec_table:: variable name ' vlist{k,1} ' is not in q2avec.name or q2avec.gname'])
     end
     if size(vlist,2)<=2
         q2avec_table(k).frcst_name = vlist{k,2};
         q2avec_table(k).plot = q2avec(indx).plot;
     else
         q2avec_table(k).frcst_name = [vlist{k,2} vlist{k,3} vlist{k,4}];
         if isequal(vlist{k,3},' Deviation')
             q2avec_table(k).plot = 2;
         else
             % isequal(vlist{k,3},' Inflation') || isequal(vlist{k,3},' Growth')
             q2avec_table(k).plot = 1;
         end
     end
     for j=1:NbFields
         q2avec_table(k).(myFields{j}) =   q2avec(indx).(myFields{j});
     end
 end
 
