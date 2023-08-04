clear all

% Read Patient from excel
data = readmatrix('DOE.xlsx','Sheet','GT','Range','D4:N62');
data(:,5:10)=[];
% L1  -  D1 - V_avg  - D0 -  FFR

% Add noise to patient data
noise = 10/100;
P_data = data(:,1:4);
FFR_GT = data(:,end);
no_p = length(P_data);

% Take average measurements to get close to true values
iter = 100;
FFR_t = zeros(no_p,1);
res_FFR_t = zeros(no_p,1);
FFR_r_t = zeros(no_p,iter);
res_FFR_r_t = zeros(no_p,iter);

for i=1:no_p
    P_noisy = Noisy_data(P_data(i,:),noise);
    FFR_r_t(i,r) = LaminarSeperation(P_noisy(1)*1e-3,P_noisy(2)*1e-3,P_noisy(3),P_noisy(4)*1e-3);   
    FFR_t(i) = LaminarSeperation(P_data(i,1)*1e-3,P_data(i,2)*1e-3,P_data(i,3),P_data(i,4)*1e-3);
    res_FFR_t(i) = FFR_t(i)-FFR_GT(i);
    for r=1:iter
        P_noisy = Noisy_data(P_data(i,:),noise);
        FFR_r_t(i,r) = LaminarSeperation(P_noisy(1)*1e-3,P_noisy(2)*1e-3,P_noisy(3),P_noisy(4)*1e-3);
        res_FFR_r_t(i,r) = FFR_r_t(i,r)-FFR_GT(i);
    end  
    
end

p_no = 10;
scatter([1:1:iter],res_FFR_r_t(p_no,:),7,'filled');
hold on
plot(res_FFR_t(p_no)*ones(1,iter))

var_t = var(res_FFR_t)
var_r_t = var(mean(res_FFR_r_t'))

function [P] = Noisy_data(P_data,noise)
P = zeros(1,4);
sigma_1 = noise*P_data(1);
sigma_2 = noise*P_data(2);
sigma_3 = noise*P_data(3);
sigma_4 = noise*P_data(4);
P(1) = P_data(1) + randn(1).*sigma_1;
P(2) = P_data(2) + randn(1).*sigma_2;
P(3) = P_data(3) + randn(1).*sigma_3;
P(4) = P_data(4) + randn(1).*sigma_4;
end