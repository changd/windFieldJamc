clear
opt = 'GFS'; 
path = '/Users/changd/Dropbox (MIT)/2019_JAMC/MATLAB_revision/';

AIC = zeros(6,5); aic2 = zeros(6,5);
logSSE = zeros(6,5); logL = zeros(6,5);
SSE = zeros(6,5); varRes = zeros(6,5);

%%
input = 'MF';
load([path, 'Results/ResultsCNLS/Results_' opt '_' input '.mat']); k = 4;
for i = 1:6
    [logSSE(1,i), AIC(1,i), SSE(1,i), varRes(1,i)] = calcAIC(models(i).res, k);
end

%%
% load([path, 'Results/ResultsCNLS/ResultsBenchmark' '.mat']); k = 6;
% for i = 1:6
%     [logSSE(2,i), AIC(2,i), SSE(2,i), varRes(2,i)] = calcAIC(models(i).resTrain, k);
% end

%%
input = 'Vst'; specs = 'p1'; 
load([path, 'Results/ResultsVst/Results_' opt '_' input '_' specs '.mat']);
k = 8;
for i = 1:6
    [logSSE(3,i), AIC(3,i), SSE(3,i), varRes(3,i)] = calcAIC(models(i).res, k);
end

%%
% input = 'Vst';  
% load([path, 'Results/ResultsVst/Results_' opt '_' input '_' 'VecAdd' '.mat']);
% k = 10;
% for i = 1:6
%     [logSSE(4,i), AIC(4,i), SSE(4,i), varRes(4,i)] = calcAIC(models(i).res, k);
% end

%%
input = 'Vsh';  
load([path, 'Results/ResultsVsh_WVN1/Results_' opt '_' input '_' 'VecAdd' '.mat']);
k = 10;
for i = 1:6
    [logSSE(5,i), AIC(5,i), SSE(5,i), varRes(5,i)] = calcAIC(models(i).res, k);
end

%%
input = 'Vst_Vsh';  
load([path, 'Results/ResultsVsh_WVN1/Results_' opt '_' input '_' 'VecAdd' '.mat']);
k = 14;
for i = 1:6
    [logSSE(6,i), AIC(6,i), SSE(6,i), varRes(6,i)] = calcAIC(models(i).res, k);
end

%%
% minAIC = min(AIC,[],2); deltaAIC = zeros(6,5);
minAIC = min(min(AIC));
for i = 1:6
    for j = 1:6
        deltaAIC(i,j) = AIC(i,j) - minAIC;
    end
end

%%
relProb = exp(-deltaAIC/2);