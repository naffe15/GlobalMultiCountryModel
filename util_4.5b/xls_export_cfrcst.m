
function xls_export_cfrcst(path, asimul, country, sheet_name)

if nargin<4
    sheet_name = 'baseline'; 
end
years = asimul(1).TimeLineQA(end-2:end);
varsT = vars_scenarios;

varsT1 = varsT;

for i = 1:length(varsT)
     if  isempty(strfind(char(varsT(i).name),'_')>0)
%     if contains(varsT(i).name,'_')!=0
varsT1(i).name = strcat(varsT(i).name, '_', country);
               
    else
varsT1(i).name= varsT(i).name; 
        
%         varsT1(i).quarterly = eval(varsT1(i).nameq);
%         varsT1(i).annual = quarterly2annual(varsT1(i).quarterly,varsT1(i).yss,varsT1(i).GYTREND0,varsT1(i).type,varsT1(i).islog,varsT1(i).aux);

    end
end



for i = 1:length(varsT1)
    %     for j = 1: length(oo_.jrc.forecastA)
%     indx = find(contains(cellstr(char(oo_.jrc.forecastA.VarName)),char(varsT1(i).name))==1);
    indx = find(strncmpi(cellstr(char(asimul.VarName)),char(varsT1(i).name),length(char(varsT1(i).name)))==1);
    if isempty(indx)
        table.(char(varsT1(i).name)).name= char(varsT1(i).name) ;
        table.(char(varsT1(i).name)).def= char(varsT1(i).def) ;
        table.(char(varsT1(i).name)).frc=  ones(3,1)*NaN;
    else
        %         if strcmp(char(varsT1(i).name),oo_.jrc.forecastA(j).VarName)==1
        table.(char(varsT1(i).name)).name= char(varsT1(i).name) ;
        table.(char(varsT1(i).name)).def= char(varsT1(i).def) ;
        table.(char(varsT1(i).name)).frc= asimul(indx(1)).exogassmpt(end-2:end);
    end
    %         end
    %     end
end

box_headers1 = {'Technical assumptions'; table.GEA_EA.name; table.INOMA_EA.name;  table.PHIBRENTOBSA_RoW.name;table.GYOBSA_RoW.name};
box_headersA1 = {'Technical assumptions'; table.GEA_EA.def; table.INOMA_EA.def;table.PHIBRENTOBSA_RoW.def; table.GYOBSA_RoW.def};
box_frc1 = [ [NaN, NaN, NaN] ;table.GEA_EA.frc'; table.INOMA_EA.frc'; table.PHIBRENTOBSA_RoW.frc';table.GYOBSA_RoW.frc'];


box_headers2 = {'Expenditure side - annual percentage changes'; table.(['GCA_' country]).name;table.(['GCGA_' country]).name;table.(['GIGA_' country]).name;table.(['GIA_' country]).name;table.(['GXA_' country]).name;table.(['GMTOTA_' country]).name;table.(['GYOBSA_' country]).name};
box_headersA2 = {'Expenditure side - annual percentage changes'; table.(['GCA_' country]).def;table.(['GCGA_' country]).def;table.(['GIGA_' country]).def;table.(['GIA_' country]).def;table.(['GXA_' country]).def;table.(['GMTOTA_' country]).def;table.(['GYOBSA_' country]).def};
box_frc2 = [[NaN, NaN, NaN] ;table.(['GCA_' country]).frc'; table.(['GCGA_' country]).frc';table.(['GIGA_' country]).frc';table.(['GIA_' country]).frc';table.(['GXA_' country]).frc';table.(['GMTOTA_' country]).frc';table.(['GYOBSA_' country]).frc'];

box_headers3 = {'Expenditure side deflators - annual percentage changes'; table.(['PHICVATA_' country]).name; table.(['PHIGA_' country]).name; table.(['PHIIGA_' country]).name; table.(['PHIIA_' country]).name; table.(['PHIXA_' country]).name; table.(['PHIMTOTA_' country]).name; table.(['PHIYOBSA_' country]).name};
box_headersA3 = {'Expenditure side deflators - annual percentage changes'; table.(['PHICVATA_' country]).def; table.(['PHIGA_' country]).def; table.(['PHIIGA_' country]).def; table.(['PHIIA_' country]).def; table.(['PHIXA_' country]).def; table.(['PHIMTOTA_' country]).def; table.(['PHIYOBSA_' country]).def};
box_frc3 = [[NaN, NaN, NaN] ;table.(['PHICVATA_' country]).frc'; table.(['PHIGA_' country]).frc'; table.(['PHIIGA_' country]).frc'; table.(['PHIIA_' country]).frc'; table.(['PHIXA_' country]).frc'; table.(['PHIMTOTA_' country]).frc'; table.(['PHIYOBSA_' country]).frc'];

box_headers4 = {'Other price indicators'; table.(['PHIWA_' country]).name};
box_headersA4 = {'Other price indicators'; table.(['PHIWA_' country]).def};
box_frc4 = [ [NaN, NaN, NaN] ;table.(['PHIWA_' country]).frc'];

box_headers5 = {'Supply of goods and services'; table.(['CUA_' country]).name; table.(['FNtNA_' country]).name};
box_headersA5 = {'Other price indicators'; table.(['CUA_' country]).def;table.(['FNtNA_' country]).def};
box_frc5 = [ [NaN, NaN, NaN] ;table.(['CUA_' country]).frc';table.(['FNtNA_' country]).frc'];

box_headers6 = {'Labour market indicators'; table.(['GNA_' country]).name};
box_headersA6 = {'Labour market indicators'; table.(['GNA_' country]).def};
box_frc6 = [[NaN, NaN, NaN] ;table.(['GNA_' country]).frc'];

% box_headers7 = {'Households+NPISH'; table.GYOBSA.name; table.INOMA.name}
% box_values7 = [[NaN, NaN, NaN] ;table.GYOBSA.value'; table.INOMA.value']

box_headers8 = {'External accounts - % of GDP'; table.(['TBYA_' country]).name};
box_headersA8 = {'External accounts - % of GDP'; table.(['TBYA_' country]).def};
box_frc8 = [[NaN, NaN, NaN] ;table.(['TBYA_' country]).frc'];

box_headers9 = {'General Government account - % of GDP'; table.(['RGYA_' country]).name;  table.(['BGYA_' country]).name;  table.(['DFETOTA_' country]).name;  table.(['GTA_' country]).name};
box_headersA9 = {'General Government account - % of GDP'; table.(['RGYA_' country]).def;table.(['BGYA_' country]).def;table.(['DFETOTA_' country]).def;table.(['GTA_' country]).def};
box_frc9 = [[NaN, NaN, NaN] ;table.(['RGYA_' country]).frc';table.(['BGYA_' country]).frc';table.(['DFETOTA_' country]).frc';table.(['GTA_' country]).frc'];

box_frcst = [box_frc1;box_frc2;box_frc3;box_frc4;box_frc5;box_frc6;box_frc8;box_frc9];

% Table preparation
table_length = length(box_headers1)+length(box_headers2)+length(box_headers3)+length(box_headers4)+length(box_headers5)+length(box_headers6)+length(box_headers8)+length(box_headers9);
box = cell(length(table_length), 10);
box(1,1) = {'Macroeconomic scenario'};
box(2,1) = {'AF18'};
box(2,2) = {'Baseline'};
box(3,1:4) = [[NaN] num2cell(years') ];
box(4:table_length+3, 1) = [box_headersA1;box_headersA2;box_headersA3;box_headersA4;box_headersA5;box_headersA6;box_headersA8;box_headersA9];
box(4:table_length+3, 6) = [box_headers1;box_headers2;box_headers3;box_headers4;box_headers5;box_headers6;box_headers8;box_headers9];
box(4:table_length+3, 2:4) = num2cell(100*box_frcst);

% writing the XLS
% sheet = [assmpt '_Cond(No assumption)'];
xlswrite(sprintf('%s\\Conditional_forecast.xls', path), box, sheet_name)

