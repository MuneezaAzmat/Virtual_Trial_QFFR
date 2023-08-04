clear all

noise = 5/100;
% Normal blood flow rate CFR=4 and normal rest bf =120ml/min , Vavg = 0.161
V_avg = 0.1592;
% V_avg = 0.161;
cutoff = 0.8 ;                                  % FFR value below which stenosis becomes significant

% Read Patient from excel
P_data = readmatrix('DOE.xlsx','Sheet','GT','Range','D4:N63');

n = 10000;
FFR_Qp = zeros(n,1);
FFR_Qn = zeros(n,1);
FFR_D  = zeros(n,1);
FFR_GT = zeros(n,1);
m = zeros(n,1);

% Loop through patients 
for subj=1:n
    prob = rand;
    if prob <0.80
        m(subj) = randi([38 60]);
    else
        m(subj) = randi([1 37]);
    end
   
    FFR_GT(subj) = P_data(m(subj),10); 
    
    % Add noise to patient data
    P_data_noisy = Noisy_data(P_data(m(subj),1:3),noise);
    D_0 = 4.6 + randn(1).*(noise*4.6);
    
    FFR_Qp(subj) = LaminarSeperation(P_data_noisy(1)*1e-3,P_data_noisy(2)*1e-3,P_data_noisy(3),D_0*1e-3);
    FFR_Qn(subj) = LaminarSeperation(P_data_noisy(1)*1e-3,P_data_noisy(2)*1e-3,V_avg,D_0*1e-3);
    FFR_D(subj)  = (P_data_noisy(2)/D_0); 
    
end
% Correction for -ive FFR
FFR_GT((FFR_GT<0))=0;
FFR_Qp((FFR_Qp<0))=0;
FFR_Qn((FFR_Qn<0))=0;

% Classification Based on FFR
Stn_GT = double(FFR_GT < cutoff);         % 1: high-risk stenosis , 0: low-risk stn 
Stn_Qp = double(FFR_Qp < cutoff);        
Stn_Qn = double(FFR_Qn < cutoff);
Stn_D  = double(FFR_D  < 0.375);

% Accuracy 
Acc_Qp = sum(double(Stn_Qp==Stn_GT))/n
Acc_Qn = sum(double(Stn_Qn==Stn_GT))/n
Acc_D  = sum(double(Stn_D ==Stn_GT))/n

% RMS error 
RMS_Qp = norm(FFR_GT-FFR_Qp)*(1/sqrt(n))
RMS_Qn = norm(FFR_GT-FFR_Qn)*(1/sqrt(n))
RMS_D  = norm(FFR_GT-FFR_D)*(1/sqrt(n))

Conf_Qp = confusionmat(Stn_GT,Stn_Qp)
Conf_Qn = confusionmat(Stn_GT,Stn_Qn)
Conf_D  = confusionmat(Stn_GT,Stn_D )

% ROC curve 
[fpr_Qp,tpr_Qp,T_Qp,auc_Qp] = perfcurve(Stn_GT,FFR_Qp,0);  %label 0 because score > 0.9 corresponds to 0
[fpr_Qn,tpr_Qn,T_Qn,auc_Qn] = perfcurve(Stn_GT,FFR_Qn,0);
[fpr_D,tpr_D,T_D,auc_D] = perfcurve(Stn_GT,FFR_D,0);

figure('visible', 'on', 'Position',[2587.4,189.8, 390.4 ,270.4])
hold on
plot(fpr_D,tpr_D, 'Color', "#77AC30",'Linewidth',2,'DisplayName',strcat('FFR_G: AUC = ',num2str(auc_D)))
plot(fpr_Qn,tpr_Qn, 'Color', "#0072BD",'Linewidth',2,'DisplayName',strcat('FFR_N: AUC = ',num2str(auc_Qn))')
plot(fpr_Qp,tpr_Qp,'r','Linewidth',2,'DisplayName',strcat('FFR_P: AUC = ',num2str(auc_Qp)))
xlabel('False Positive rate','fontsize',14)
ylabel('True Positive rate','fontsize',14)
%title ('ROC curves for classification of stenosis based on FFR')
legend ('fontsize',12,'Location','southeast')
grid on

function [P] = Noisy_data(P_data,noise)
    P = zeros(1,3);
    sigma_1 = noise*P_data(1);
    sigma_2 = noise*P_data(2);
    sigma_3 = noise*P_data(3);
    P(1) = P_data(1) + randn(1).*sigma_1;
    P(2) = P_data(2) + randn(1).*sigma_2;
    P(3) = P_data(3) + randn(1).*sigma_3;
end
