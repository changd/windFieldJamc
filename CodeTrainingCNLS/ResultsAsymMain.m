%% MF
clear
opt = 'GFS';
disp('MF')
noWvn = 1;
x0 = zeros(8,1);
b_l = zeros(8,1);
b_u = zeros(8,1);

tic
[x, SSE_train, SSE_val, SSE_IreneHW, SSE_KatrinaHW, SSE_SandyHW, SSE_IsabelHW, ...
SSE_KatrinaWRF, SSE_SandyWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt);
toc

path = ['/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/'];
save([path, 'Results/ResultsMF.mat'],'x','SSE_train','SSE_val');
save([path, 'Results/ResultsTestMF.mat'],'SSE_IreneHW','SSE_KatrinaHW','SSE_SandyHW', ...
    'SSE_IsabelHW');
% save('Results/ResultsTestWRF_MF.mat','SSE_KatrinaWRF','SSE_SandyWRF');

%% Benchmark
clear
opt = 'GFS';

[SSE_train, SSE_val, SSE_IreneHW, SSE_KatrinaHW, SSE_SandyHW, ...
        SSE_IsabelHW, SSE_KatrinaWRF, SSE_SandyWRF] = ResultsBenchmark(opt);
    
path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
save([path, 'Results/ResultsBenchmark.mat'],'SSE_train','SSE_val');
save([path, 'Results/ResultsTestBenchmark.mat'],'SSE_IreneHW','SSE_KatrinaHW','SSE_SandyHW', ...
    'SSE_IsabelHW');   
% save('Results/ResultsTestWRF_MF.mat','SSE_KatrinaWRF','SSE_SandyWRF');

%% CNLS Vst p=0
clear
opt = 'GFS'; method = 'CNLS'; input = 'Vst'; specs = 'p0'; noWvn = 1;
disp('CNLS Vst p=0');
x0 = [2.93; 0.39; 89.5*(pi/180); 3.4*(pi/180); zeros(4,1)];
b_l = [0;0;-10;0; zeros(4,1)]; b_u = [10;0;10;0; zeros(4,1)];

tic
[x, SSE_train, SSE_val, SSE_IreneHW, SSE_KatrinaHW, SSE_SandyHW, SSE_IsabelHW, ...
SSE_KatrinaWRF, SSE_SandyWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt);
toc

path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
save([path, 'Results/Results_' opt '_' method '_' input '_' specs '.mat'],'x','SSE_train','SSE_val');
save([path, 'Results/ResultsTest_' opt '_' method '_' input '_' specs '.mat'],'SSE_IreneHW','SSE_KatrinaHW','SSE_SandyHW', ...
    'SSE_IsabelHW');
% save(['Results/ResultsTestWRF_' opt '_' method '_' input '_' specs '.mat'],'SSE_KatrinaWRF','SSE_SandyWRF');

%% CNLS Vst p=1
clear
opt = 'GFS'; method = 'CNLS'; input = 'Vst'; specs = 'p1'; noWvn = 1;
disp('CNLS Vst p=1')
bl_St = [0;-10;-10;-10]; bu_St = [10;10;10;10];
x0 = [2.93; 0.39; 89.5*(pi/180); 3.4*(pi/180); zeros(4,1)];
b_l=[bl_St; zeros(4,1)]; b_u=[bu_St; zeros(4,1)];

tic
[x, SSE_train, SSE_val, SSE_IreneHW, SSE_KatrinaHW, SSE_SandyHW, SSE_IsabelHW, ...
SSE_KatrinaWRF, SSE_SandyWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt);
toc

path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
save([path, 'Results/Results_' opt '_' method '_' input '_' specs '.mat'],'x','SSE_train','SSE_val');
save([path, 'Results/ResultsTest_' opt '_' method '_' input '_' specs '.mat'],'SSE_IreneHW','SSE_KatrinaHW','SSE_SandyHW', ...
    'SSE_IsabelHW');
% save(['Results/ResultsTestWRF_' opt '_' method '_' input '_' specs '.mat'],'SSE_KatrinaWRF','SSE_SandyWRF');

%%
clear
opt = 'GFS';
xSt = [2.93; 0.39; 89.5*(pi/180); 3.4*(pi/180)];
xShr1 = zeros(4,1);
xShr2 = [1.41; 0.03; -60.2*(pi/180); 1.8*(pi/180)];
xShr3 = [0.93; -0.00; -25.1*(pi/180); 0.1*(pi/180)];

bl_St = [0;-10;-10;-10]; bl_Shr0 = [0;0;-10;0]; bl_Shr1 = [0;-10;-10;-10];

bu_St = [10;10;10;10]; bu_Shr0 = [10;0;10;0]; bu_Shr1 = [10;10;10;10];

%% CNLS Vsh Wvn-2
opt = 'GFS'; method = 'CNLS'; input = 'Vsh'; specs = 'Wvn2';
disp('CNLS Sh WVN-2')

inputInit = 'Vst'; specsInit = 'p1';
path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
load([path, 'Results/Results_' opt '_' method '_' inputInit '_' specsInit '.mat']);
noWvn = 2;
x0St = x(1:4,:); x0Sh = repmat([xShr1; xShr2],1,5);
x0=[xSt; zeros(4,1); 0.1; 0; -45*(pi/180); 0];
b_l=[bl_St; zeros(4,1); repmat(bl_Shr0,[noWvn-1,1])];
b_u=[bu_St; zeros(4,1); repmat(bu_Shr0,[noWvn-1,1])];

%%
tic
[x, SSE_train, SSE_val, SSE_IreneHW, SSE_KatrinaHW, SSE_SandyHW, SSE_IsabelHW, ...
SSE_KatrinaWRF, SSE_SandyWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt);
toc

path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
save([path, 'Results/Results_' opt '_' method '_' input '_' specs '.mat'],'x','SSE_train','SSE_val');
save([path, 'Results/ResultsTest_' opt '_' method '_' input '_' specs '.mat'],'SSE_IreneHW','SSE_KatrinaHW','SSE_SandyHW', ...
    'SSE_IsabelHW');
% save(['Results/ResultsTestWRF_' opt '_' method '_' input '_' specs '.mat'],'SSE_KatrinaWRF','SSE_SandyWRF');

%% CNLS Vsh Wvn-2+3
opt = 'GFS'; method = 'CNLS'; input = 'Vsh'; specs = 'Wvn3';
disp('CNLS Sh WVN-2+3')

inputInit = 'Vsh'; specsInit = 'Wvn2';
path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
load([path, 'Results/Results_' opt '_' method '_' inputInit '_' specsInit '.mat']);
noWvn = 3;
xSt = mean(x,2);
% x0=[xSt(1:8); 0.1; 0; -45*(pi/180); 0; 0.1; 0; -45*(pi/180); 0];
x0 = [mean(x,2); 0.1; 0; -45*(pi/180); 0];
b_l=[bl_St; zeros(4,1); repmat(bl_Shr0,[noWvn-1,1])];
b_u=[bu_St; zeros(4,1); repmat(bu_Shr0,[noWvn-1,1])];

%%
tic
[x, SSE_train, SSE_val, SSE_IreneHW, SSE_KatrinaHW, SSE_SandyHW, SSE_IsabelHW, ...
SSE_KatrinaWRF, SSE_SandyWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt);
toc

path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
save([path, 'Results/Results_' opt '_' method '_' input '_' specs '.mat'],'x','SSE_train','SSE_val');
save([path, 'Results/ResultsTest_' opt '_' method '_' input '_' specs '.mat'],'SSE_IreneHW','SSE_KatrinaHW','SSE_SandyHW', ...
    'SSE_IsabelHW');
% save(['Results/ResultsTestWRF_' opt '_' method '_' input '_' specs '.mat'],'SSE_KatrinaWRF','SSE_SandyWRF');