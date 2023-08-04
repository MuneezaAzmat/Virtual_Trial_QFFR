clear all 

% Read Patient from excel
P_data = readmatrix('DOE.xlsx','Sheet','GT','Range','D4:M62');
P_data(:,5:9)=[];
% L1  -  D1 - V_avg  - D0 - %Dstn - FFR 

% Test-Train split 30-70%
%test_id = randi([1 59],1,17);
test_id = [ 55   4  38   2   1  29  32   4  35  10   7  24  57  23  30  10  19]';
train_id = setdiff([1:1:59]', test_id);

% Add noise to patient data
N = 10000;
Train_data  = zeros(N,8); %6 features from excel and one label
cutoff=0.8;

for subj=1:N
    % Generate data with noise and prevalence
    prob = rand;
    if prob <0.80
        m(subj) = train_id(randi([25 44]));
    else                                           
        m(subj) = train_id(randi([1 24]));
    end
    
    noise = randi([5 10])/100;
    Train_data(subj,1:3) = Noisy_data(P_data(m(subj),1:3),noise);
    D_0 = 4.6 + randn(1).*(noise*4.6);
    Train_data(subj,4) =  D_0;
    Train_data(subj,5) =  (1-(Train_data(subj,2)/D_0))*100;
    Train_data(subj,6) = LaminarSeperation(Train_data(subj,1)*1e-3,Train_data(subj,2)*1e-3,Train_data(subj,3),Train_data(subj,4)*1e-3);
    Train_data(subj,7) = P_data(m(subj),end); % true FFR
    Train_data(subj,8) = double(Train_data(subj,7) < cutoff); % true label
end

%% Generate Test data

n = 2500;
Test_data  = zeros(n,8);
for subj=1:n
    id = test_id(randi([1 17]));
    noise = randi([5 15])/100;
    
    Test_data(subj,1:3) = Noisy_data(P_data(m(subj),1:3),noise);
    Test_data(subj,4) =  4.6 + randn(1).*(noise*4.6);
    Test_data(subj,5) =  (1-(Test_data(subj,2)/D_0))*100;
    Test_data(subj,6) = LaminarSeperation(Test_data(subj,1)*1e-3,Test_data(subj,2)*1e-3,Test_data(subj,3),Test_data(subj,4)*1e-3);
    Test_data(subj,7) = P_data(m(subj),end); % true FFR
    Test_data(subj,8) = double(Test_data(subj,7) < cutoff); % true label
end


function [P] = Noisy_data(P_data,noise)
    P = zeros(1,3);
    
    sigma_1 = noise*P_data(1);
    sigma_2 = noise*P_data(2);
    sigma_3 = noise*P_data(3);
    P(1) = P_data(1) + randn(1).*sigma_1;
    P(2) = P_data(2) + randn(1).*sigma_2;
    P(3) = P_data(3) + randn(1).*sigma_3;
end