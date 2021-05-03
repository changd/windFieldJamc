%%
clear
opt = 'GFS'; input = 'Vsh';
path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';

%% Entire data set (WVN-1+2)
input = 'Vsh'; specs = 'WVN1_2'; noWvn = 2; disp('CNLS Vsh WVN-1+2')
Rbound = [0 Inf];
x0 = repmat([2.93; 0.39; 89.5*(pi/180); 3.4*(pi/180); 3.66; 0.02; -31.9*(pi/180); ...
    -2.8*(pi/180); 0.1;0;-45*(pi/180);0], [1 5]);
b_l = [0;-10;-10;-10; 0;-10;-10;-10; 0;-10;-10;-10];
b_u = [10;10;10;10; 10;10;10;10; 10;10;10;10];

[models, MSE_TestHWind, MSE_TestWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound);
save([path, 'Results/ResultsVsh_WVN2_3/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'MSE_TestHWind', 'MSE_TestWRF');

%% Entire data set (WVN-2)
input = 'Vsh'; specs = 'WVN2'; noWvn = 2; disp('CNLS Vsh WVN-2')
Rbound = [0 Inf];
x0 = repmat([2.93; 0.39; 89.5*(pi/180); 3.4*(pi/180); 3.66; 0.02; -31.9*(pi/180); ...
    -2.8*(pi/180); 0.1;0;-45*(pi/180);0], [1 5]);
b_l = [0;-10;-10;-10; 0;0;0;0; 0;-10;-10;-10];
b_u = [10;10;10;10; 0;0;0;0; 10;10;10;10];

[models, MSE_TestHWind, MSE_TestWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound);
save([path, 'Results/ResultsVsh_WVN2_3/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'MSE_TestHWind', 'MSE_TestWRF');

%% Entire data set (WVN-2+3)
input = 'Vsh'; specs = 'WVN1_2_3'; noWvn = 3; disp('CNLS Vsh WVN-2+3')
Rbound = [0 Inf];
x0 = repmat([2.93; 0.39; 89.5*(pi/180); 3.4*(pi/180); 3.66; 0.02; -31.9*(pi/180); ...
    -2.8*(pi/180); 0.1;0;-45*(pi/180);0; 0.1;0;-45*(pi/180);0], [1 5]);
b_l = [0;-10;-10;-10; 0;0;0;0; 0;-10;-10;-10; 0;-10;-10;-10];
b_u = [10;10;10;10; 0;0;0;0; 10;10;10;10; 10;10;10;10];

[models, MSE_TestHWind, MSE_TestWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound);
save([path, 'Results/ResultsVsh_WVN2_3/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'MSE_TestHWind', 'MSE_TestWRF');

%% Entire data set (WVN-3)
input = 'Vsh'; specs = 'WVN3'; noWvn = 3; disp('CNLS Vsh WVN-3')
Rbound = [0 Inf];
x0 = repmat([2.93; 0.39; 89.5*(pi/180); 3.4*(pi/180); 3.66; 0.02; -31.9*(pi/180); ...
    -2.8*(pi/180); 0.1;0;-45*(pi/180);0; 0.1;0;-45*(pi/180);0], [1 5]);
b_l = [0;-10;-10;-10; 0;0;0;0; 0;0;0;0; 0;-10;-10;-10];
b_u = [10;10;10;10; 0;0;0;0; 0;0;0;0; 10;10;10;10];

[models, MSE_TestHWind, MSE_TestWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound);
save([path, 'Results/ResultsVsh_WVN2_3/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'MSE_TestHWind', 'MSE_TestWRF');

%%
specsInit = 'WVN1';
load([path, 'Results/ResultsVsh_WVN1/Results_' opt '_' input '_' specsInit '.mat']);

xInit = zeros(8,5);
for i = 1:length(models)
    xInit(:,i) = models(i).x;
end

%% Entire data set (WVN-1+2+3)
input = 'Vsh'; specs = 'WVN1_2_3'; noWvn = 3; disp('CNLS Vsh WVN-1+2+3')
Rbound = [0 Inf];
x0 = [xInit; repmat([0.1;0;-45*(pi/180);0; 0.1;0;-45*(pi/180);0], [1 5])];
b_l = [0;-10;-10;-10; 0;-10;-10;-10; 0;-10;-10;-10; 0;-10;-10;-10];
b_u = [10;10;10;10; 10;10;10;10; 10;10;10;10; 10;10;10;10];

[models, MSE_TestHWind, MSE_TestWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound);
save([path, 'Results/ResultsVsh_WVN2_3/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'MSE_TestHWind', 'MSE_TestWRF');

%% Entire data set (WVN-1+3)
input = 'Vsh'; specs = 'WVN1_3'; noWvn = 3; disp('CNLS Vsh WVN-1+3')
Rbound = [0 Inf];
x0 = [xInit; repmat([0.1;0;-45*(pi/180);0; 0.1;0;-45*(pi/180);0], [1 5])];
b_l = [0;-10;-10;-10; 0;-10;-10;-10; 0;0;0;0; 0;-10;-10;-10];
b_u = [10;10;10;10; 10;10;10;10; 0;0;0;0; 10;10;10;10];

[models, MSE_TestHWind, MSE_TestWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound);
save([path, 'Results/ResultsVsh_WVN2_3/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'MSE_TestHWind', 'MSE_TestWRF');