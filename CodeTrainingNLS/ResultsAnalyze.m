clear
opt = 'GFS'; method = 'CNLS'; input = 'Vst'; specs = 'p1'; noWvn = 1;
path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
load([path, 'Results/ResultsNLS/Results_' opt '_' method '_' input '_' specs '.mat']);

%%
i = 1;
R = models(i).R_Train;
res = models(i).resTrain; resSq = res.^2;
% idx = find(R >= 1);
idx = find(resSq < 15);
resFilter = res(idx);
resSqFilter = resFilter.^2; 
hist(res);