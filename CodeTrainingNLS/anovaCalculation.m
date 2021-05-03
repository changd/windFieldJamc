clear
opt = 'GFS'; 
path = '/Users/changd/Dropbox (MIT)/2019_JAMC/MATLAB_revision/';

AIC = zeros(6,5); aic2 = zeros(6,5);
logSSE = zeros(6,5); logL = zeros(6,5);
SSE = zeros(6,5); varRes = zeros(6,5);
resStruct = struct([]);

%%
input = 'MF';
load([path, 'Results/ResultsCNLS/Results_' opt '_' input '.mat']); k = 4;
for i = 1:length(models)
    resStruct(i).resMF = models(i).res;
end

%%
load([path, 'Results/ResultsCNLS/ResultsBenchmark' '.mat']); k = 6;
for i = 1:length(models)
    resStruct(i).resTV = models(i).resTrain;
end

%%
input = 'Vst'; specs = 'p1'; 
load([path, 'Results/ResultsVst/Results_' opt '_' input '_' specs '.mat']); k = 4;
for i = 1:length(models)
    resStruct(i).resWVN1 = models(i).res;
end

%%
input = 'Vst';  
load([path, 'Results/ResultsVst/Results_' opt '_' input '_' 'VecAdd' '.mat']);
for i = 1:length(models)
    resStruct(i).resTV_WVN1 = models(i).res;
end

%%
input = 'Vsh';  
load([path, 'Results/ResultsVsh_WVN1/Results_' opt '_' input '_' 'VecAdd' '.mat']);
for i = 1:length(models)
    resStruct(i).resShr = models(i).res;
end

%%
input = 'Vst_Vsh';  
load([path, 'Results/ResultsVsh_WVN1/Results_' opt '_' input '_' 'VecAdd' '.mat']);
for i = 1:length(models)
    resStruct(i).resWVN1_Shr = models(i).res;
end

%%
for i = 1:length(models)
    resStruct(i).anova_p(1) = anova1([resStruct(i).resMF resStruct(i).resTV], [], 'off');
    resStruct(i).anova_p(2) = anova1([resStruct(i).resTV resStruct(i).resWVN1], [], 'off');
    resStruct(i).anova_p(3) = anova1([resStruct(i).resWVN1 resStruct(i).resTV_WVN1], [], 'off');
    resStruct(i).anova_p(4) = anova1([resStruct(i).resTV resStruct(i).resTV_WVN1], [], 'off');
    resStruct(i).anova_p(5) = anova1([resStruct(i).resTV_WVN1 resStruct(i).resWVN1_Shr], [], 'off');
    resStruct(i).anova_p(6) = anova1([resStruct(i).resTV resStruct(i).resShr], [], 'off');
    resStruct(i).anova_p(7) = anova1([resStruct(i).resWVN1_Shr resStruct(i).resShr], [], 'off');
end