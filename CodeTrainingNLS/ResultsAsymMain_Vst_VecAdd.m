clear
opt = 'GFS'; 
path = '/Users/changd/Dropbox (MIT)/2019_JAMC/MATLAB_revision/';

%%
input = 'Vst'; noWvn = 1; disp('CNLS Vst p=1')
bl_St = [-5;-10;-10;-10]; bu_St = [10;10;10;10]; Rbound = [0 Inf];
x0 = repmat([4.56; 0.07; -100*(pi/180); 4.52*(pi/180); 3.37; 0.13; -11.78*(pi/180); ...
    -4.34*(pi/180)],[1,6]);
b_l=[bl_St; zeros(4,1)]; b_u=[bu_St; zeros(4,1)];
% b_l(3) = 80*pi/180; b_u(3) = 100*pi/180;

models = ResultsTrainTestVecAdd(x0, b_l, b_u, noWvn, opt, Rbound);

%%
save([path, 'Results/ResultsVst/Results_' opt '_' input '_VecAdd' '.mat'], 'models', 'MSE_TestWRF');

%% Radial bins
Rbounds = [0.8 1.2; 1.8 2.2; 2.9 3.2; 3.9 4.2; 4.9 5.2];
[nBounds, ~] = size(Rbounds);

for i = 1:nBounds
    disp(i)
    input = 'Vst'; specs = ['p1_Rm_' num2str(i)]; noWvn = 1; disp('CNLS Vst p=1 (Rm)')
    bl_St = [-5;-10;-10;-10]; bu_St = [10;10;10;10]; Rbound = Rbounds(i,:);
    x0 = repmat([4.56; 0.07; 9.85*(pi/180); 4.52*(pi/180); 3.37; 0.13; -11.78*(pi/180); ...
    -4.34*(pi/180)],[1,5]);
    b_l=[bl_St; zeros(4,1)]; b_u=[bu_St; zeros(4,1)];

    [models, SSE_TestHWind, SSE_TestWRF] = ResultsTrainTestVecAdd(x0, b_l, b_u, noWvn, opt, Rbound);
    save([path, 'Results/ResultsVst/Results_' opt '_' input '_VecAdd' '.mat'],'models', 'SSE_TestHWind', 'SSE_TestWRF');
end

%% Storm-by-storm
opt = 'SHIPS_Rev'; 
x0 = [4.56; 0.07; 9.85*(pi/180); 4.52*(pi/180); 3.37; 0.13; -11.78*(pi/180); ...
    -4.34*(pi/180)];
b_l=[bl_St; zeros(4,1)]; b_u=[bu_St; zeros(4,1)]; 
noWvn = 1; Rbounds = [0 Inf];
b_l(1) = 0; 

% Dennis
stormName = 'Dennis'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVst/Results_' stormName '_' opt '_' input '_VecAdd' '.mat'], 'x', 'ci', 'MSE');

% Ingrid
stormName = 'Ingrid'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVst/Results_' stormName '_' opt '_' input '_VecAdd' '.mat'], 'x', 'ci', 'MSE');

% Isaac
stormName = 'Isaac'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVst/Results_' stormName '_' opt '_' input '_VecAdd' '.mat'], 'x', 'ci', 'MSE');

% Sandy
stormName = 'SandyHW_Train_'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVst/Results_' stormName '_' opt '_' input '_VecAdd' '.mat'], 'x', 'ci', 'MSE');

% Wilma
stormName = 'Wilma'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVst/Results_' stormName '_' opt '_' input '_VecAdd' '.mat'], 'x', 'ci', 'MSE');