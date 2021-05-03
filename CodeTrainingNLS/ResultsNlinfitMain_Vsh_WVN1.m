% function ResultsNlinfitMain(opt, path)
clear
opt = 'GFS'; 
path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';

%% CNLS Vst p=1
input = 'Vst'; specs = 'p1'; noWvn = 1; modelType = 'p0';
disp('CNLS Vst p=1');

filePathInitial = [path, 'Results/ResultsCNLS/Results_' opt '_' input '_' specs '.mat'];
load(filePathInitial); beta0 = x(1:4,:);

[models, beta] = ResultsNlinfit(beta0, 'GFS', modelType);
save([path, 'Results/ResultsNLS/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'beta');

%% CNLS Vsh Wvn-1 p=1 (Use CAE initialization for Vtr and Vsh)
opt = 'GFS'; input = 'Vsh'; specs = 'Wvn1_p1_VtrCAE_VshCAE'; modelType = 'p1';
disp('CNLS Sh WVN-1 p=1')

filePathInitial = [path, 'Results/ResultsCNLS/Results_' opt '_' input '_' specs '.mat'];
load(filePathInitial);
beta0 = x; %beta0 = [beta0(1:4,:); beta0(9:12,:)];

[models, beta] = ResultsNlinfit(beta0, 'GFS', modelType);
save([path, 'Results/ResultsNLS/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'beta');

%% CNLS Vsh Wvn-1 p=1 (Use CAE initialization for Vtr, U14 for Vsh)
opt = 'GFS'; input = 'Vsh'; specs = 'Wvn1_p1_VtrCAE_VshU14'; modelType = 'p1';
disp('CNLS Sh WVN-1 p=1')

filePathInitial = [path, 'Results/ResultsCNLS/Results_' opt '_' input '_' specs '.mat'];
load(filePathInitial);
beta0 = x; %beta0 = [beta0(1:4,:); beta0(9:12,:)];

[models, beta] = ResultsNlinfit(beta0, 'GFS', modelType);
save([path, 'Results/ResultsNLS/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'beta');

%% CNLS Vsh Wvn-1 p=1 (Use U14 initialization for Vtr and Vsh)
opt = 'GFS'; input = 'Vsh'; specs = 'Wvn1_p1_VtrU14_VshU14'; modelType = 'p1';
disp('CNLS Sh WVN-1 p=1')

filePathInitial = [path, 'Results/ResultsCNLS/Results_' opt '_' input '_' specs '.mat'];
load(filePathInitial);
beta0 = x; %beta0 = [beta0(1:4,:); beta0(9:12,:)];

%%
[models, beta] = ResultsNlinfit(beta0, 'GFS', modelType);
save([path, 'Results/ResultsNLS/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'beta');

%% CNLS Vsh Wvn-1 p=1 (Use U14 initialization for Vtr, CAE for Vsh)
opt = 'GFS'; input = 'Vsh'; specs = 'Wvn1_p1_VtrU14_VshCAE'; modelType = 'p1';
disp('CNLS Sh WVN-1 p=1')

filePathInitial = [path, 'Results/ResultsCNLS/Results_' opt '_' input '_' specs '.mat'];
load(filePathInitial);
beta0 = x; %beta0 = [beta0(1:4,:); beta0(9:12,:)];

[models, beta] = ResultsNlinfit(beta0, 'GFS', modelType);
save([path, 'Results/ResultsNLS/Results_' opt '_' input '_' specs '.mat'], ...
    'models', 'beta');