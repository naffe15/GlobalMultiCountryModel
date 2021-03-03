function dataobs_change(fname,oo_,varobs,nfrcst,exo_path,occbin)

% global oo_ 

if nargin<5 || isempty(exo_path) || isempty(occbin)
    cases=0;
end
if nargin<6 || isempty(exo_path)
    cases = 0;
end
if nargin<7 || isempty(occbin)
    cases = 0;
end

if exo_path ==1, cases = 1; end
if occbin ==1, cases = 2; end

%saveas(fname,[])
copyfile([fname '.mat'],[ fname '_old.mat']);

eval(['load ' fname])
vv = who;
% ww = fieldnames(oo_.forecast.Mean); 
ww = varobs; 
index = find(ismember(vv,ww));

string = ['save ',fname,' T '];
for jj = 1 : length(index)
    if sum(isnan(eval(vv{index(jj)})))<10 % there are not too many nans else estimation check fails
        if cases == 0
            eval([	vv{index(jj)} '(end-nfrcst+1:end) = NaN;'  ])    
        elseif cases ==1
            eval([	vv{index(jj)} '(end-nfrcst+1:end) = oo_.forecast_exo_path.Mean.' vv{index(jj)} ';'  ])    
        elseif cases ==2 
            eval([	vv{index(jj)} '(end-nfrcst+1:end) = oo_.occbin_forecast.Mean.' vv{index(jj)} ';'  ])    
        end
    end
    string = [ string ' ' vv{index(jj)}];    
end

eval(string);

