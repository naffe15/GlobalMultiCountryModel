
function compare_scenarios2(path, base, s1, country, assmpt)

if nargin<5
    assmpt = 'new';
end
years = base.oo_.jrc.forecastA(1).TimeLineQA(end-2:end);
varsT = vars_scenarios2;
% country = 'EA';
% assmpt = 'INOM_EA_High';
% % varsT = struct('name',{'INOMA', 'GYOBSA'}, 'def',{'interest', 'GDP'})
% years = [2017 2018 2019];

varsT1 = varsT;
% load(q2avec);
% idir = findstr('_',assmpt);
% temp = eraseBetween(assmpt,idir(end),length(assmpt));
% assmpt1 = q2avec(find(contains(char(q2avec.qname),temp)==1)).name;

for i = 1:length(varsT)
%     Cambiare la istruzione per gli ZEPS
%        if ismember( contains(varsT(i).name,'ZEPS'), contains(varsT(i).name,'REA') )==0
        
        if (  (contains(varsT(i).name,'ZEPS' ) )== 1 && ~( (contains(varsT(i).name,'REA' ) ) || (contains(varsT(i).name,'RoW' ) ))  )
            varsT1(i).name = strcat(varsT(i).name, '_', country); 
   
     else 
         
      if ((contains(varsT(i).name,'ZEPS' ) )== 1 && ( (contains(varsT(i).name,'REA' ) ) || (contains(varsT(i).name,'RoW' ) )) )

      varsT1(i).name= varsT(i).name; 
   
      else 
         
      if ((contains(varsT(i).name,'ZEPS' ) )== 1 && ( (contains(varsT(i).name,'IT_IT' ) ) )) 

      varsT1(i).name= varsT(i).name; 
   
      else  
          
         
        if  isempty(strfind(char(varsT(i).name),'_')>0)
         
%     if contains(varsT(i).name,'_')!=0
        varsT1(i).name = strcat(varsT(i).name, '_', country);
               
        else
        
        varsT1(i).name= varsT(i).name; 
        
%         varsT1(i).quarterly = eval(varsT1(i).nameq);
%         varsT1(i).annual = quarterly2annual(varsT1(i).quarterly,varsT1(i).yss,varsT1(i).GYTREND0,varsT1(i).type,varsT1(i).islog,varsT1(i).aux);
        end
         end
    end
     end
end



for i = 1:length(varsT1)
    %     for j = 1: length(oo_.jrc.forecastA)
%     indx = find(contains(cellstr(char(oo_.jrc.forecastA.VarName)),char(varsT1(i).name))==1);
    indx = find(strncmpi(cellstr(char(base.oo_.jrc.forecastA.VarName)),char(varsT1(i).name),length(char(varsT1(i).name)))==1);
    if isempty(indx)
        table.(char(varsT1(i).name)).value= ones(3,1)*NaN ;
        table.(char(varsT1(i).name)).name= char(varsT1(i).name) ;
        table.(char(varsT1(i).name)).def= char(varsT1(i).def) ;
        table.(char(varsT1(i).name)).frc=  ones(3,1)*NaN;
    else
        %         if strcmp(char(varsT1(i).name),oo_.jrc.forecastA(j).VarName)==1
        table.(char(varsT1(i).name)).value= s1.oo_.jrc.forecastA(indx(1)).exogassmpt(end-2:end) ;
        table.(char(varsT1(i).name)).name= char(varsT1(i).name) ;
        table.(char(varsT1(i).name)).def= char(varsT1(i).def) ;
        table.(char(varsT1(i).name)).frc= base.oo_.jrc.forecastA(indx(1)).exogassmpt(end-2:end);
    end
    %         end
    %     end
end

box_headers1 = {'TECHNICAL ASSUMPTIONS'; table.GEA_EA.name; table.INOMA_EA.name;  table.PHIBRENTOBSA_RoW.name;table.GYOBSA_RoW.name};
box_headersA1 = {'TECHNICAL ASSUMPTIONS'; table.GEA_EA.def; table.INOMA_EA.def;table.PHIBRENTOBSA_RoW.def; table.GYOBSA_RoW.def};
box_values1 = [ [NaN, NaN, NaN] ;table.GEA_EA.value'; table.INOMA_EA.value';table.PHIBRENTOBSA_RoW.value'; table.GYOBSA_RoW.value'];
box_frc1 = [ [NaN, NaN, NaN] ;table.GEA_EA.frc'; table.INOMA_EA.frc'; table.PHIBRENTOBSA_RoW.frc';table.GYOBSA_RoW.frc'];

box_headers2 = {'DEMAND SIDE - ANNUAL PERCENTAGE CHANGES'; ...
    table.(['GYOBSA_' country]).name;...
    table.(['GCA_' country]).name;...
    table.(['GCSA_' country]).name;...
    table.(['GCCA_' country]).name;...
    table.(['GCGA_' country]).name;...
    table.(['GTA_' country]).name;...
    table.(['GIGA_' country]).name;...
    table.(['GIA_' country]).name;...
    table.(['GXA_' country]).name;...
    table.(['GMTOTA_' country]).name};

box_headersA2 = {'DEMAND SIDE - ANNUAL PERCENTAGE CHANGES'; ...
    table.(['GYOBSA_' country]).def;...
    table.(['GCA_' country]).def;...
    table.(['GCSA_' country]).def;...
    table.(['GCCA_' country]).def;...
    table.(['GCGA_' country]).def;...
    table.(['GTA_' country]).def;...
    table.(['GIGA_' country]).def;...
    table.(['GIA_' country]).def;...
    table.(['GXA_' country]).def;...
    table.(['GMTOTA_' country]).def};

box_values2 = [[NaN, NaN, NaN] ; ...
    table.(['GYOBSA_' country]).value';...
    table.(['GCA_' country]).value';...
    table.(['GCSA_' country]).value';...
    table.(['GCCA_' country]).value';...
    table.(['GCGA_' country]).value';...
    table.(['GTA_' country]).value';...
    table.(['GIGA_' country]).value';...
    table.(['GIA_' country]).value';...
    table.(['GXA_' country]).value';...
    table.(['GMTOTA_' country]).value'];

box_frc2 = [ [NaN, NaN, NaN] ; ...
     table.(['GYOBSA_' country]).frc';...
    table.(['GCA_' country]).frc'; ...
    table.(['GCSA_' country]).frc';...
    table.(['GCCA_' country]).frc';...
    table.(['GCGA_' country]).frc';...
    table.(['GTA_' country]).frc';...
    table.(['GIGA_' country]).frc';...
    table.(['GIA_' country]).frc';...
    table.(['GXA_' country]).frc';...
    table.(['GMTOTA_' country]).frc'];


box_headers3 = {'DEMAND DEFLATORS - ANNUAL PERCENTAGE CHANGES';...
    table.(['PHIYOBSA_' country]).name;...
    table.(['PHICVATA_' country]).name; ...
    table.(['PHIGA_' country]).name; ...
    table.(['PHIIGA_' country]).name; ...
    table.(['PHIIA_' country]).name; ...
    table.(['PHIXA_' country]).name;...
    table.(['PHIMTOTA_' country]).name};

box_headersA3 = {'DEMAND DEFLATORS - ANNUAL PERCENTAGE CHANGES'; ...
    table.(['PHIYOBSA_' country]).def;...
    table.(['PHICVATA_' country]).def; ...
    table.(['PHIGA_' country]).def; ...
    table.(['PHIIGA_' country]).def; ...
    table.(['PHIIA_' country]).def; ...
    table.(['PHIXA_' country]).def;...
    table.(['PHIMTOTA_' country]).def};

box_values3 = [[NaN, NaN, NaN] ;...
    table.(['PHIYOBSA_' country]).value';...
    table.(['PHICVATA_' country]).value';...
    table.(['PHIGA_' country]).value'; ...
    table.(['PHIIGA_' country]).value';...
    table.(['PHIIA_' country]).value';...
    table.(['PHIXA_' country]).value';...
    table.(['PHIMTOTA_' country]).value'];

box_frc3 = [[NaN, NaN, NaN] ;...
    table.(['PHIYOBSA_' country]).frc';...
    table.(['PHICVATA_' country]).frc';...
    table.(['PHIGA_' country]).frc';...
    table.(['PHIIGA_' country]).frc';...
    table.(['PHIIA_' country]).frc'; ...
    table.(['PHIXA_' country]).frc';...
    table.(['PHIMTOTA_' country]).frc'];

% % %% 
% % % % % % % % % % % OTHER PRICE INDICATORS
% % % % % % 
    box_headers4 = {'OTHER PRICE INDICATORS'; ...

     table.(['PHIWA_' country]).name; ...
     table.(['INOMA_' country]).name;...
     table.(['INOMGA_' country]).name; ...
    table.(['INOMGLTA_' country]).name; ...
    table.(['RA_' country]).name; ...
     table.(['RSA_' country]).name; ...
     table.(['GPYMRKUPA_' country]).name; ...
    table.(['INOMSA_' country]).name};

    box_headersA4 = {'OTHER PRICE INDICATORS'; ...
        
     table.(['PHIWA_' country]).def; ...
     table.(['INOMA_' country]).def;...
    table.(['INOMGA_' country]).def;...
    table.(['INOMGLTA_' country]).def;...
    table.(['RA_' country]).def;...
    table.(['RSA_' country]).def;...
    table.(['GPYMRKUPA_' country]).def;...
    table.(['INOMSA_' country]).def};


    box_values4 = [ [NaN, NaN, NaN] ; ...
        
     table.(['PHIWA_' country]).value'; ...
     table.(['INOMA_' country]).value';...
     table.(['INOMGA_' country]).value';...
     table.(['INOMGLTA_' country]).value';...
    table.(['RA_' country]).value';...
   table.(['RSA_' country]).value';...
   table.(['GPYMRKUPA_' country]).value';...
    table.(['INOMSA_' country]).value'] ; 

    box_frc4 = [ [NaN, NaN, NaN] ; ...
        
     table.(['PHIWA_' country]).frc'; ...
    table.(['INOMA_' country]).frc';...
    table.(['INOMGA_' country]).frc';...
    table.(['INOMGLTA_' country]).frc';...
    table.(['RA_' country]).frc';...
    table.(['RSA_' country]).frc';...
    table.(['GPYMRKUPA_' country]).frc';...
    table.(['INOMSA_' country]).frc'] ; 

box_headers5 = {'SUPPLY OF GOODS AND SERVICES'; table.(['CUA_' country]).name; table.(['FNtNA_' country]).name};
box_headersA5 = {'SUPPLY OF GOODS AND SERVICES'; table.(['CUA_' country]).def;table.(['FNtNA_' country]).def};
box_values5 = [ [NaN, NaN, NaN] ;table.(['CUA_' country]).value';table.(['FNtNA_' country]).value'];
box_frc5 = [ [NaN, NaN, NaN] ;table.(['CUA_' country]).frc';table.(['FNtNA_' country]).frc'];

    box_headers6 = {'LABOUR MARKET INDICATORS';...
    table.(['GNA_' country]).name;...
    table.(['YLA_' country]).name};
    
    box_headersA6 = {'LABOUR MARKET INDICATORS';...
    table.(['GNA_' country]).def;...
    table.(['YLA_' country]).def};

    box_values6 = [[NaN, NaN, NaN] ;...
    table.(['GNA_' country]).value';...
    table.(['YLA_' country]).value'];
  
    box_frc6 = [[NaN, NaN, NaN] ;.....
    table.(['GNA_' country]).frc';...
    table.(['YLA_' country]).frc'];

% box_headers7 = {'Households+NPISH'; table.GYOBSA.name; table.INOMA.name}
% box_values7 = [[NaN, NaN, NaN] ;table.GYOBSA.value'; table.INOMA.value']

    box_headers8 = {'EXTERNAL ACCOUNTS ';...
    table.(['NFAYA_' country]).name;...
    table.(['TBYA_' country]).name};

    box_headersA8 = {'EXTERNAL ACCOUNTS';...
    table.(['NFAYA_' country]).def;...    
    table.(['TBYA_' country]).def};

     box_values8 = [[NaN, NaN, NaN] ;....
     table.(['NFAYA_' country]).value' ;...   
     table.(['TBYA_' country]).value'];
    
     box_frc8 = [[NaN, NaN, NaN] ;...
    table.(['NFAYA_' country]).frc';...
    table.(['TBYA_' country]).frc'];



box_headers9 = {'GENERAL GOVERNMENT ACCOUNT '; ...
    table.(['BGYA_' country]).name; ...
    table.(['PSGYA_' country]).name ; ...
    table.(['IBYA_' country]).name; ...
    table.(['RGYA_' country]).name; ...
    table.(['RGTAUNYA_' country]).name; ...
    table.(['RGTAUCYA_' country]).name; ...
    table.(['RGTAUKYA_' country]).name; ...
    table.(['RGTAXYA_' country]).name; ... 
    table.(['RGTAUOILYA_' country]).name;...
    table.(['DFETOTA_' country]).name;...
    table.(['CGYA_' country]).name;...
    table.(['TYA_' country]).name}; 

    box_headersA9 = {'GENERAL GOVERMENT ACCOUNT'; ...
    table.(['BGYA_' country]).def; ...
    table.(['PSGYA_' country]).def ; ...
    table.(['IBYA_' country]).def ; ...
    table.(['RGYA_' country]).def; ...
    table.(['RGTAUNYA_' country]).def ;...
    table.(['RGTAUCYA_' country]).def;...
    table.(['RGTAUKYA_' country]).def;...
    table.(['RGTAXYA_' country]).def;...
    table.(['RGTAUOILYA_' country]).def;...
    table.(['DFETOTA_' country]).def; ...
    table.(['CGYA_' country]).def; ...
    table.(['TYA_' country]).def };  

    box_values9 = [[NaN, NaN, NaN] ;...
    table.(['BGYA_' country]).value';...
    table.(['PSGYA_' country]).value';...
    table.(['IBYA_' country]).value' ; ...
    table.(['RGYA_' country]).value'; ...
    table.(['RGTAUNYA_' country]).value';...
    table.(['RGTAUCYA_' country]).value'; ...
    table.(['RGTAUKYA_' country]).value'; ...
    table.(['RGTAXYA_' country]).value';...
    table.(['RGTAUOILYA_' country]).value'; ...
    table.(['DFETOTA_' country]).value'; ...
    table.(['CGYA_' country]).value'; ...
    table.(['TYA_' country]).value' ];  

    box_frc9 = [[NaN, NaN, NaN] ;...
    table.(['BGYA_' country]).frc'; 
    table.(['PSGYA_' country]).frc';...
    table.(['IBYA_' country]).frc';...
    table.(['RGYA_' country]).frc';...
    table.(['RGTAUNYA_' country]).frc';...
    table.(['RGTAUCYA_' country]).frc';...
    table.(['RGTAUKYA_' country]).frc'; 
    table.(['RGTAXYA_' country]).frc';...
    table.(['RGTAUOILYA_' country]).frc';...
    table.(['DFETOTA_' country]).frc';...
    table.(['CGYA_' country]).frc';...
    table.(['TYA_' country]).frc'];


     box_headers10 = {'CHECK on SHOCKS - DOMESTIC ';...
     table.(['ZEPS_APC_' country]).name;...
     table.(['ZEPS_APG_' country]).name;...
     table.(['ZEPS_API_' country]).name;...
     table.(['ZEPS_APIG_' country]).name;...
     table.(['ZEPS_AY_' country]).name;...
     table.(['ZEPS_CU_' country]).name;...
     table.(['ZEPS_FQ_' country]).name;...
     table.(['ZEPS_FN_' country]).name;...
     table.(['ZEPS_ND_' country]).name;...
     table.(['ZEPS_BW_' country]).name;...
     table.(['ZEPS_B_' country]).name;...
     table.(['ZEPS_G_' country]).name;...
     table.(['ZEPS_ACTR_' country]).name;...
     table.(['ZEPS_HPERE_' country]).name;...
     table.(['ZEPS_PARTR_' country]).name;...
     table.(['ZEPS_POP_' country]).name;...
     table.(['ZEPS_IG_' country]).name;...
      
     table.(['ZEPS_PX_' country]).name;...
     table.(['ZEPS_S_' country]).name;...
     table.(['ZEPS_T_' country]).name;...
     table.(['ZEPS_TAUC_' country]).name;...
     table.(['ZEPS_TC_' country]).name;...
     table.(['ZEPS_TAUN_' country]).name;...
     table.(['ZEPS_TAUK_' country]).name;...
     table.(['ZEPS_TAX_' country]).name;...
     table.(['ZEPS_U_' country]).name;...
     table.(['ZEPS_UC_' country]).name;...
     table.(['ZEPS_BEN_' country]).name};
      
 
     box_headersA10 = {'CHECK on SHOCKS - DOMESTIC '; ...
      table.(['ZEPS_APC_' country]).def;... 
      table.(['ZEPS_APG_' country]).def;...
      table.(['ZEPS_API_' country]).def;...
      table.(['ZEPS_APIG_' country]).def;...
      table.(['ZEPS_AY_' country]).def;...
      table.(['ZEPS_CU_' country]).def;...
      table.(['ZEPS_FQ_' country]).def;...
      table.(['ZEPS_FN_' country]).def;...
      table.(['ZEPS_ND_' country]).def;...
      table.(['ZEPS_BW_' country]).def;...
      table.(['ZEPS_B_' country]).def;...
      table.(['ZEPS_G_' country]).def;...
      table.(['ZEPS_ACTR_' country]).def;...
      table.(['ZEPS_HPERE_' country]).def;...
      table.(['ZEPS_PARTR_' country]).def;...
      table.(['ZEPS_POP_' country]).def;...
      table.(['ZEPS_IG_' country]).def;...
      
      table.(['ZEPS_PX_' country]).def;...
      table.(['ZEPS_S_' country]).def;...
      table.(['ZEPS_T_' country]).def;...
      table.(['ZEPS_TAUC_' country]).def;...
      table.(['ZEPS_TC_' country]).def;...
      table.(['ZEPS_TAUN_' country]).def;...
      table.(['ZEPS_TAUK_' country]).def;...
      table.(['ZEPS_TAX_' country]).def;...
      table.(['ZEPS_U_' country]).def;...
      table.(['ZEPS_UC_' country]).def;...
      table.(['ZEPS_BEN_' country]).def};
     
 
     box_values10 = [[NaN, NaN, NaN] ;...
     table.(['ZEPS_APC_' country]).value'; ...   
     table.(['ZEPS_APG_' country]).value'; ... 
     table.(['ZEPS_API_' country]).value'; ...
     table.(['ZEPS_APIG_' country]).value'; ...
     table.(['ZEPS_AY_' country]).value'; ...
     table.(['ZEPS_CU_' country]).value'; ...
     table.(['ZEPS_FQ_' country]).value';...
      table.(['ZEPS_FN_' country]).value';...
      table.(['ZEPS_ND_' country]).value';...
      table.(['ZEPS_BW_' country]).value';...
      table.(['ZEPS_B_' country]).value';...
      table.(['ZEPS_G_' country]).value';...
      table.(['ZEPS_ACTR_' country]).value';...
      table.(['ZEPS_HPERE_' country]).value';...
      table.(['ZEPS_PARTR_' country]).value';...
      table.(['ZEPS_POP_' country]).value';...
      table.(['ZEPS_IG_' country]).value';...
      
      table.(['ZEPS_PX_' country]).value';...
      table.(['ZEPS_S_' country]).value';...
      table.(['ZEPS_T_' country]).value';...
      table.(['ZEPS_TAUC_' country]).value';...
      table.(['ZEPS_TC_' country]).value';...
      table.(['ZEPS_TAUN_' country]).value';...
      table.(['ZEPS_TAUK_' country]).value';...
      table.(['ZEPS_TAX_' country]).value';...
      table.(['ZEPS_U_' country]).value';...
      table.(['ZEPS_UC_' country]).value';...
      table.(['ZEPS_BEN_' country]).value'];
      
 
     box_frc10 = [[NaN, NaN, NaN] ;...
     table.(['ZEPS_APC_' country]).frc';...  
     table.(['ZEPS_APG_' country]).frc';...
     table.(['ZEPS_API_' country]).frc';...
     table.(['ZEPS_APIG_' country]).frc';...
     table.(['ZEPS_AY_' country]).frc';...
     table.(['ZEPS_CU_' country]).frc';...
     table.(['ZEPS_FQ_' country]).frc';...
     table.(['ZEPS_FN_' country]).frc';...
     table.(['ZEPS_ND_' country]).frc';...
     table.(['ZEPS_BW_' country]).frc';...
     table.(['ZEPS_B_' country]).frc';...
     table.(['ZEPS_G_' country]).frc';...
     table.(['ZEPS_ACTR_' country]).frc';...
     table.(['ZEPS_HPERE_' country]).frc';...
     table.(['ZEPS_PARTR_' country]).frc';...
     table.(['ZEPS_POP_' country]).frc';...
     table.(['ZEPS_IG_' country]).frc';...
      
      table.(['ZEPS_PX_' country]).frc';...
      table.(['ZEPS_S_' country]).frc';...
      table.(['ZEPS_T_' country]).frc';...
      table.(['ZEPS_TAUC_' country]).frc';...
      table.(['ZEPS_TC_' country]).frc';...
      table.(['ZEPS_TAUN_' country]).frc';...
      table.(['ZEPS_TAUK_' country]).frc';...
      table.(['ZEPS_TAX_' country]).frc';...
      table.(['ZEPS_U_' country]).frc';...
      table.(['ZEPS_UC_' country]).frc';...
      table.(['ZEPS_BEN_' country]).frc'];
     
  
  
  %%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&&&&&&&&&&&&&REA
  
  
   box_headers11 = {' REST OF THE EA (REA) BLOCK ';...
     
   
     table.(['GYOBSA_REA']).name;...
     table.(['GCA_REA']).name;...
     table.(['PHIYOBSA_REA']).name;...
     table.(['INOMA_REA']).name;...
     table.(['TBYA_REA']).name;...  
   
     table.(['ZEPS_POP_REA']).name;...
     table.(['ZEPS_UC_REA']).name;...
     table.(['ZEPS_M_REA' ]).name;...
     table.(['ZEPS_Y_REA' ]).name;...
     table.(['ZEPS_PX_REA']).name;...
     table.(['ZEPS_A_REA']).name...
      };
      
 
     box_headersA11 = {'REST OF THE EA (REA) BLOCK  '; ...
      
      
     
     table.(['GYOBSA_REA']).def;...
     table.(['GCA_REA']).def;...
     table.(['PHIYOBSA_REA']).def;...
     table.(['INOMA_REA']).def;...
     table.(['TBYA_REA']).def;...  
     
     
      table.(['ZEPS_POP_REA']).def;...
      table.(['ZEPS_UC_REA']).def;...
      table.(['ZEPS_M_REA']).def;...
      table.(['ZEPS_Y_REA']).def;...
      table.(['ZEPS_PX_REA']).def;...
      table.(['ZEPS_A_REA']).def...
      
      };
     
 
     box_values11 = [[NaN, NaN, NaN] ;...
    
      
     table.(['GYOBSA_REA']).value';...
     table.(['GCA_REA']).value';...
     table.(['PHIYOBSA_REA']).value';...
     table.(['INOMA_REA']).value';...
     table.(['TBYA_REA']).value';...
     
     
      table.(['ZEPS_POP_REA' ]).value';...
      table.(['ZEPS_UC_REA']).value';...
      table.(['ZEPS_M_REA']).value';...
      table.(['ZEPS_Y_REA']).value';...
      table.(['ZEPS_PX_REA']).value';...
      table.(['ZEPS_A_REA']).value'...
      
      ];
      
 
     box_frc11 = [[NaN, NaN, NaN] ;...
     
      
     table.(['GYOBSA_REA']).frc';...
     table.(['GCA_REA']).frc';...
     table.(['PHIYOBSA_REA']).frc';...
     table.(['INOMA_REA']).frc';...
     table.(['TBYA_REA']).frc';...
     
      table.(['ZEPS_POP_REA' ]).frc';...
      table.(['ZEPS_UC_REA']).frc';...
      table.(['ZEPS_M_REA']).frc';...
      table.(['ZEPS_Y_REA']).frc';...
      table.(['ZEPS_PX_REA']).frc';...
      table.(['ZEPS_A_REA']).frc'...
      
      ];
     
    

     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ROW
  
  
    box_headers12 = {' REST of THE WORLD (RoW) BLOCK  ';...
        
    table.(['GYOBSA_RoW']).name;...
     table.(['GCA_RoW']).name;...
     table.(['PHIYOBSA_RoW']).name;...
     table.(['INOMA_RoW']).name;...
     table.(['TBYA_RoW']).name;...  
  
     table.(['ZEPS_POP_RoW']).name;...
     table.(['ZEPS_UC_RoW']).name;...
     table.(['ZEPS_M_RoW' ]).name;...
     table.(['ZEPS_Y_RoW' ]).name;...
     table.(['ZEPS_INOM_RoW']).name;...
     table.(['ZEPS_A_RoW']).name...
     };
      
 
     box_headersA12 = {'REST of THE WORLD (RoW) BLOCK '; ...
      
     table.(['GYOBSA_RoW']).def;...
     table.(['GCA_RoW']).def;...
     table.(['PHIYOBSA_RoW']).def;...
     table.(['INOMA_RoW']).def;...
     table.(['TBYA_RoW']).def;...  
      
      table.(['ZEPS_POP_RoW']).def;...
      table.(['ZEPS_UC_RoW']).def;...
      table.(['ZEPS_M_RoW']).def;...
      table.(['ZEPS_Y_RoW']).def;...
      table.(['ZEPS_INOM_RoW']).def;...
      table.(['ZEPS_A_RoW']).def...
      };
     
 
     box_values12 = [[NaN, NaN, NaN] ;...
    
     
      table.(['GYOBSA_RoW']).value';...
      table.(['GCA_RoW']).value';...
     table.(['PHIYOBSA_RoW']).value';...
     table.(['INOMA_RoW']).value';...
     table.(['TBYA_RoW']).value';...
     
      
      table.(['ZEPS_POP_RoW' ]).value';...
      table.(['ZEPS_UC_RoW']).value';...
      table.(['ZEPS_M_RoW']).value';...
      table.(['ZEPS_Y_RoW']).value';...
      table.(['ZEPS_INOM_RoW']).value';...
      table.(['ZEPS_A_RoW']).value'...
      ];
      
 
     box_frc12 = [[NaN, NaN, NaN] ;...
     
     
     table.(['GYOBSA_RoW']).frc';...
     table.(['GCA_RoW']).frc';...
     table.(['PHIYOBSA_RoW']).frc';...
     table.(['INOMA_RoW']).frc';...
     table.(['TBYA_RoW']).frc';...
     
      table.(['ZEPS_POP_RoW' ]).frc';...
      table.(['ZEPS_UC_RoW']).frc';...
      table.(['ZEPS_M_RoW']).frc';...
      table.(['ZEPS_Y_RoW']).frc';...
      table.(['ZEPS_INOM_RoW']).frc';...
      table.(['ZEPS_A_RoW']).frc'...
      ];
  
  
 
 %%%%%%%%%PUT ALL TOGETHER
 box_frcst = [box_frc1;box_frc2;box_frc3;box_frc4;box_frc9;box_frc10;box_frc11;box_frc12];

 box_frcst = [box_frc1;box_frc2;box_frc3;box_frc4;box_frc5;box_frc6;box_frc8;box_frc9;box_frc10;box_frc11;box_frc12];

% box_assmpt = [box_values1;box_values2;box_values3;box_values4;box_values9];

 box_assmpt = [box_values1;box_values2;box_values3;box_values4;box_values5;box_values6;box_values8;box_values9;box_values10;box_values11;box_values12];

box_dif = box_assmpt-box_frcst;

% Table preparation
% table_length = length(box_headers1)+length(box_headers2)+length(box_headers3)+length(box_headers4)+length(box_headers9);

 table_length = length(box_headers1)+length(box_headers2)+length(box_headers3)+length(box_headers4)+length(box_headers5)+length(box_headers6)+length(box_headers8)+length(box_headers9)+length(box_headers10)+length(box_headers11)+length(box_headers12);
box = cell(length(table_length), 10);
box(1,1) = {'Macroeconomic scenario'};
box(2,1) = {'AF18'};
box(2,2) = {'Baseline'};
box(2,5) = {assmpt};
box(2,8) = {'Differences'};
box(3,1:10) = [[NaN] num2cell(years') num2cell(years') num2cell(years')];
% Setting the order of the above blocks to be displayed
% box(4:table_length+3, 1) = [box_headersA1;box_headersA2;box_headersA3;box_headersA4;box_headersA9];

 box(4:table_length+3, 1) = [box_headersA1;box_headersA2;box_headersA3;box_headersA4;box_headersA5;box_headersA6;box_headersA8;box_headersA9;box_headersA10;box_headersA11;box_headersA12];

% box(4:table_length+3, 11) = [box_headers1;box_headers2;box_headers3;box_headers4;box_headers9];
 box(4:table_length+3, 11) = [box_headers1;box_headers2;box_headers3;box_headers4;box_headers5;box_headers6;box_headers8;box_headers9;box_headers10;box_headers11;box_headers12];

box(4:table_length+3, 2:4) = num2cell(100*box_frcst);
box(4:table_length+3, 5:7) = num2cell(100*box_assmpt);
box(4:table_length+3, 8:10) = num2cell(100*box_dif);

% writing the XLS
% sheet = [assmpt '_Cond(No assumption)'];
if length(assmpt)>31
xlswrite(sprintf('%s\\Compare_Scenarios2.xls', path), box, assmpt(1:31))
else
xlswrite(sprintf('%s\\Compare_Scenarios2.xls', path), box, assmpt)
end

