clear all

n=2500;
FFR_Q = zeros(n,1);
FFR_G  = zeros(n,1);
load('TrainData_varnoise')
load('TestData_varnoise')

FFR_GT = Test_data(:,end);

for i=1:n
FFR_Q(i) = LaminarSeperation(Test_data(i,1)*1e-3,Test_data(i,2)*1e-3,Test_data(i,3),Test_data(i,4)*1e-3);
FFR_G(i)  = (Test_data(i,2)/Test_data(i,4));
end
% Correction for -ive FFR
FFR_Q((FFR_Q<0))=0;

RMSE_Q = sqrt(sum((FFR_Q - FFR_GT).^2)/n);
RMSE_G  = sqrt(sum((FFR_G - FFR_GT).^2)/n);


%% Train the ML models 
[trainedModelQ, RMSE_ML_Q] = RQ_GPR(Train_data);
[trainedModel, RMSE_ML_G] = RQ_GPR_Geometry(Train_data);

% remove true ffr values
Test_data(:,end)=[];
FFR_ML_Q = trainedModelQ.predictFcn(Test_data); 
% remove flow data
Test_data(:,3)=[];
FFR_ML_G = trainedModel.predictFcn(Test_data); 


% % Scatter plot 
% 
% figure(1)
% scatter(FFR_GT,FFR_D,'filled','jitter', 'on', 'jitterAmount', 0.01);
% axis('equal')
% grid on
% hold on 
% plot(FFR_GT,FFR_GT,'k')
% xlabel('True FFR value')
% ylabel('Predicted FFR value')
% title('Geometric model')
% 
% figure(2)
% scatter(FFR_GT,FFR_Qp,'filled','jitter', 'on', 'jitterAmount', 0.01);
% axis('equal')
% grid on
% hold on 
% plot(FFR_GT,FFR_GT,'k')
% xlabel('True FFR value')
% ylabel('Predicted FFR value')
% title('Flow informed analytic model')


%% Bland Altman %%
corrinfo = {'n','r2','eq','rho'}; % stats to display of correlation scatter plot
BAinfo = {}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
label = {'True FFR','Predicted FFR'};
gnames={'FFR'};
if 1 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [0 0 1;... % or RGB triplets
		      1 0 0];
end

[cr_Q, fig_Q, stats_Q] = BlandAltman(FFR_GT, FFR_Q,label,'Flow informed analytic model',gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors,'markerSize',2);
[cr_G, fig_G, stats_G] = BlandAltman(FFR_GT, FFR_G,label,'Geometry informed analytical model',gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors,'markerSize',2);
[cr_ML_Q, fig_ML_Q, stats_ML_Q] = BlandAltman(FFR_GT, FFR_ML_Q,label,'Flow informed Machine learning model',gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors,'markerSize',2);
[cr_ML_G, fig_ML_G, stats_ML_G] = BlandAltman(FFR_GT, FFR_ML_G,label,'Geometric Machine leaning model',gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors,'markerSize',2);


