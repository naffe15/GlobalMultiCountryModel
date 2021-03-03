function [best_ar,best_ma,param_ar,param_ma,se,res,irf] = best_arma(y)

% This function identifies the order of the ARMA model
% It uses as an input a stationary time series y
% and gives the orders of the process, the estimated parameters
% the standard error of residuals and the residuals

addpath z:\SVN\captain\captain\

[p_1,q_1,C_1,D_1,Pcc_1,resvar_1,eef_1,RR_1]=ivarmaid(y/std(y),[1 0 4 4],10,1);
[p_2,q_2,C_2,D_2,Pcc_2,resvar_2,eef_2,RR_2]=ivarmaid(y/std(y),[1 0 4 4],10,2);
[p_3,q_3,C_3,D_3,Pcc_3,resvar_3,eef_3,RR_3]=ivarmaid(y/std(y),[1 0 4 4],10,3);

for i=1:size(RR_1,1)
model_1(i,1) = str2num(strcat(num2str(RR_1(i,1)),num2str(RR_1(i,2))));
end

for i=1:size(RR_1,1)
model_2(i,1) = str2num(strcat(num2str(RR_2(i,1)),num2str(RR_2(i,2))));
end

for i=1:size(RR_1,1)
model_3(i,1) = str2num(strcat(num2str(RR_3(i,1)),num2str(RR_3(i,2))));
end

rank_1     = [[1:1:size(RR_1,1)]' model_1];
[~,idx]    = sort(rank_1(:,2));
rank_1_ord = rank_1(idx,:);

rank_2     = [[1:1:size(RR_2,1)]' model_2];
[~,idx]    = sort(rank_2(:,2));
rank_2_ord = rank_2(idx,:);

rank_3     = [[1:1:size(RR_3,1)]' model_3];
[~,idx]    = sort(rank_3(:,2));
rank_3_ord = rank_3(idx,:);

rank_total    = rank_1_ord(:,1)+rank_2_ord(:,1)+rank_3_ord(:,1);
rank_min      = min(rank_total);
[best_rank,~] = find(rank_total==rank_min);
best_model    = rank_1_ord(best_rank(1),2)

best_model_string = num2str(best_model);
best_ar           = str2num(best_model_string(1));
best_ma           = str2num(best_model_string(2));

    if size(best_model_string,2) == 3
           best_ar = str2num(best_model_string(1:2));
           best_ma = str2num(best_model_string(3));
    end

[best_ar,best_ma,param_ar,param_ma,Pcc_best,se,res,RR_best] = ivarmaid(y,[best_ar best_ma best_ar best_ma],10,3);

    %if max(abs(roots(param_ar))) >= 1
    %        error('The model is explosive')
    %end
    
    %if max(abs(roots(param_ma))) >= 1
    %       error('The model is not invertible')
    %end
   
    NbModels = length(rank_total);
    k=1;
    
    while(k<NbModels)
    
        if max(abs(roots(param_ar))) >= 1 | max(abs(roots(param_ma))) >= 1

            sort_rank_total   = sort(rank_total);
            [best_rank,~]     = find(rank_total==sort_rank_total(k));
            best_model        = rank_1_ord(best_rank(1),2);

            best_model_string = num2str(best_model);
            best_ar           = str2num(best_model_string(1));
            best_ma           = str2num(best_model_string(2));
        else
            k = NbModels;
        end
        k=k+1;
        
    end
    
    if size(best_model_string,2) == 3
           best_ar = str2num(best_model_string(1:2));
           best_ma = str2num(best_model_string(3));
    end

[best_ar,best_ma,param_ar,param_ma,Pcc_best,se,res,RR_best]=ivarmaid(y,[best_ar best_ma best_ar best_ma],10,3);



% if mod(best_model,10) == 0
%  
%         clearvars best_model;
%         C=aic(y,4,1);
%         best_model = num2str(strcat(num2str(length(C)-1),num2str(0)));
%         [th,res,y0]=mar(y,length(C)-1);
%         [C]=getpar(th);
%         [param_ar]=getpar(th);
%         [param_ma]=1;
%         se=std(res);
%         disp('roots of AR part')
%         abs(roots(C));
%         %if max(abs(roots(C))) >= 1
%         %    error('The model is explosive')
%         %end
%         
%         best_model_string=num2str(best_model);
%         best_ar=str2num(best_model_string(1));
%         best_ma=str2num(best_model_string(2));
%         
%         if size(best_model_string,2)==3
%             best_ar=str2num(best_model_string(1:2));
%             best_ma=str2num(best_model_string(3));
%         end
%         
% end

irf=filter(param_ma,param_ar,[1; zeros(100,1)]);

end

