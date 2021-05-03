%%
clear
opt = 'GFS'; 
path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
x0 = repmat([2.93; 0.39; 89.5*(pi/180); 3.4*(pi/180); 3.66; 0.02; -31.9*(pi/180); ...
    -2.8*(pi/180)], [1 5]);

%% Entire data set
input = 'Vsh'; specs = 'WVN1'; noWvn = 1; disp('CNLS Vsh WVN-1')
Rbound = [0 Inf];
x0 = repmat([4.56; 0.07; 9.85*(pi/180); 4.52*(pi/180); 3.37; 0.13; -11.78*(pi/180); ...
    -4.34*(pi/180)], [1 6]);
b_l = [0;-10;-10;-10; 0;-10;-10;-10];
b_u = [10;10;10;10; 10;10;10;10];

models = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound);

%%
save([path, 'Results/ResultsVsh_WVN1/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'MSE_TestHWind', 'MSE_TestWRF');

%% Radial bins
Rbounds = [0.9 1.2; 1.8 2.2; 2.8 3.2; 3.8 4.2; 4.8 5.2];
[nBounds, ~] = size(Rbounds);

for i = 1:nBounds
    input = 'Vsh'; specs = ['WVN1_Rm_' num2str(i)]; noWvn = 1; disp(['CNLS Vsh WVN-1 ' num2str(i)])
    Rbound = Rbounds(i,:);
    x0 = repmat([4.56; 0.07; 9.85*(pi/180); 4.52*(pi/180); 3.37; 0.13; -11.78*(pi/180); ...
    -4.34*(pi/180)], [1 5]);
    b_l = [0;-10;-10;-10; 0;-10;-10;-10];
    b_u = [10;10;10;10; 10;10;10;10];
    
    [models, MSE_TestHWind, MSE_TestWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound);
    save([path, 'Results/ResultsVsh_WVN1/Results_' opt '_' input '_' specs '.mat'],'models', 'MSE_TestHWind', 'MSE_TestWRF');
end

%% Storm-by-storm
x0 = [4.56; 0.07; 9.85*(pi/180); 4.52*(pi/180); 3.37; 0.13; -11.78*(pi/180); ...
    -4.34*(pi/180)];
b_l = [0;-10;-10;-10; 0;-10;-10;-10];
b_u = [10;10;10;10; 10;10;10;10];
noWvn = 1; Rbounds = [0 Inf];
opt = 'SHIPS'; 

% Dennis
stormName = 'Dennis'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVsh_WVN1/Results_' stormName '_' opt '_' input '_' specs '.mat'], 'x', 'ci', 'MSE');

% Ingrid
stormName = 'Ingrid'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVsh_WVN1/Results_' stormName '_' opt '_' input '_' specs '.mat'], 'x', 'ci', 'MSE');

% Isaac
stormName = 'Isaac'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVsh_WVN1/Results_' stormName '_' opt '_' input '_' specs '.mat'], 'x', 'ci', 'MSE');

% Sandy
stormName = 'SandyHW_Train_'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVsh_WVN1/Results_' stormName '_' opt '_' input '_' specs '.mat'], 'x', 'ci', 'MSE');

% Wilma
stormName = 'Wilma'; disp(stormName); 
load([path, 'Training Data/TrainingStorms/', stormName, opt, '.mat']);
[x, ci, MSE] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds);
save([path, 'Results/ResultsVsh_WVN1/Results_' stormName '_' opt '_' input '_' specs '.mat'], 'x', 'ci', 'MSE');
