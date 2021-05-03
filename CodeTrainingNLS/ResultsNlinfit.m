function [models, beta] = ResultsNlinfit(beta0, opt, modelType)
    
    tic
    models = getResultsTrain(beta0, opt, modelType);
    disp('Get Results')
    toc

    tic
    models = ResultsTestHWind(models, opt, modelType);
    disp('Testing HWind')
    toc
    
    beta = zeros(length(models(1).beta), 5);
    for i = 1:length(models)
        beta(:,i) = models(i).beta;
    end
%     
% %     [SSE_KatrinaWRF, SSE_SandyWRF] = ResultsTestWRF(x, noWvn, opt);
% %     disp('Testing WRF')
%     SSE_KatrinaWRF = 0; SSE_SandyWRF = 0;
end

function models = getResultsTrain(beta0, opt, modelType)
    idx = 1:1:5;
%     x = zeros(4+4*noWvn,length(idx));
%     SSE_train = zeros(1,length(idx));
%     SSE_val = zeros(1,length(idx));
    path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
    
    models = struct([]);
    for i = 1:5
%         i
        filePath = [path, 'Training Data/TrainingSets/Training', ...
            num2str(i), '_', opt, '.mat'];
        models = getResults(models, i, filePath, beta0(:,i), modelType);
        models = getResidualStatisticsTraining(models, i, filePath);
    end
end

function models = getResults(models, index, filePath, beta0, modelType)

    load(filePath);

    Storm = StormTrain;
    [X, V, radiusTrain] = getInputs(Storm);
    if strcmp(modelType, 'p0') == 1
        [beta, R, J, CovB, MSE, ~] = nlinfit(X, V, @velocityModel_p0, beta0);
    elseif strcmp(modelType, 'p1') == 1
        [beta, R, J, CovB, MSE, ~] = nlinfit(X, V, @velocityModel_p1, beta0);
    elseif strcmp(modelType, 'p01') == 1
        [beta, R, J, CovB, MSE, ~] = nlinfit(X, V, @velocityModel_p01, beta0);
    elseif strcmp(modelType, 'p10') == 1
        [beta, R, J, CovB, MSE, ~] = nlinfit(X, V, @velocityModel_p10, beta0);
    end
    CI = nlparci(beta, R, 'covar', CovB);
    
    Storm = StormVal;
    [X, V, radiusVal] = getInputs(Storm);
    
    if strcmp(modelType, 'p0') == 1
        Y = velocityModel_p0(beta, X);
    elseif strcmp(modelType, 'p1') == 1
        Y = velocityModel_p1(beta, X);
    elseif strcmp(modelType, 'p01') == 1
        Y = velocityModel_p01(beta, X);
    elseif strcmp(modelType, 'p10') == 1
        Y = velocityModel_p10(beta, X);
    end
    R_Val = V - Y;
    
    models(index).beta = beta;
    models(index).resTrain = R;
    models(index).resVal = R_Val;
    models(index).J = J;
    models(index).covB = CovB;
    models(index).MSE = MSE;
    models(index).CI = CI;
    models(index).R_Train = radiusTrain;
    models(index).R_Val = radiusVal;
end

function models = ResultsTestHWind(models, opt, modelType)

    path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
    
    for i = 1:5
        load([path 'Training Data/TestStorms/IreneHW_' opt '.mat']);
        [models(i).resIrene, models(i).R_Irene] = ResultsTest(Storm, models(i).beta, modelType);
        [models(i).MSE_Irene, models(i).prctile16_Irene, models(i).prctile50_Irene, ...
        models(i).prctile84_Irene, models(i).resRadialBin_Irene] = ...
        getResidualStatistics(models(i).resIrene, Storm);
        
        load([path 'Training Data/TestStorms/KatrinaHW_' opt '.mat']);
        [models(i).resKatrina, models(i).R_Katrina] = ResultsTest(Storm, models(i).beta, modelType);
        [models(i).MSE_Katrina, models(i).prctile16_Katrina, models(i).prctile50_Katrina, ...
        models(i).prctile84_Katrina, models(i).resRadialBin_Katrina] = ...
        getResidualStatistics(models(i).resKatrina, Storm);
        
        load([path 'Training Data/TestStorms/SandyHW_' opt '.mat']);
        [models(i).resSandy, models(i).R_Sandy] = ResultsTest(Storm, models(i).beta, modelType);
        [models(i).MSE_Sandy, models(i).prctile16_Sandy, models(i).prctile50_Sandy, ...
        models(i).prctile84_Sandy, models(i).resRadialBin_Sandy] = ...
        getResidualStatistics(models(i).resSandy, Storm);
    
        load([path 'Training Data/TestStorms/IsabelHW_SHIPS.mat']);
        [models(i).resIsabel, models(i).R_Isabel] = ResultsTest(Storm, models(i).beta, modelType);
        [models(i).MSE_Isabel, models(i).prctile16_Isabel, models(i).prctile50_Isabel, ...
        models(i).prctile84_Isabel, models(i).resRadialBin_Isabel] = ...
        getResidualStatistics(models(i).resIsabel, Storm);        
    end
end

function [Res, R] = ResultsTest(Storm, beta, modelType)
    [X, V, R] = getInputs(Storm);
    
    if strcmp(modelType, 'p0') == 1
        Y = velocityModel_p0(beta, X);
    elseif strcmp(modelType, 'p1') == 1
        Y = velocityModel_p1(beta, X);
    elseif strcmp(modelType, 'p01') == 1
        Y = velocityModel_p01(beta, X);
    elseif strcmp(modelType, 'p10') == 1
        Y = velocityModel_p10(beta, X);
    end
    
    Res = V - Y;
end

function [X, V, R] = getInputs(Storm)
    R = [];
    Y = []; Vm = []; dirTrans = []; Vtr = []; dirShr = []; Vsh = []; V = [];
    [Rcut, ~, ~, ~, ~, ~] = CutoffCalc2(Storm);
    for i = 1:length(Storm)
        [YIter, ~, idx] = VHolCNLS_Iter(Storm(i).AsymR, Storm(i).Vm, ...
            Storm(i).Rmax, Storm(i).B, Storm(i).Rn, Storm(i).Xn, Rcut(i));

        L = length(YIter);
        Y = [Y; YIter'];
        Vm = [Vm; Storm(i).Vm*ones(L,1)];
        dirTrans = [dirTrans; Storm(i).AsymDir(idx)'];
        Vtr = [Vtr; Storm(i).Vtrans*ones(L,1)];
        dirShr = [dirShr; Storm(i).AsymDirShr(idx)'];
        Vsh = [Vsh; Storm(i).Vshr*ones(L,1)];
        V = [V; Storm(i).AsymV(idx)'];
        R = [R; Storm(i).AsymR(idx)'/Storm(i).Rmax];
    end
    X = [Y Vm dirTrans Vtr dirShr Vsh];
end

function [Y, R, idx] = VHolCNLS_Iter(R, Vm, Rmax, B, Rn, Xn, Rcut)
    
    Xcalc = @(Rm, Rn, Xn, R) 0.5+(R-Rm)*(Xn-0.5)/(Rn-Rm);
    V_HRev=@(Vm,Rm,B,R,N,X) ((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^X;
    
    idx=find(R < Rcut); 
    R = R(idx);

    X = zeros(1, length(R));
    for j = 1:length(R)
        if R(j) > Rmax
            X(j) = Xcalc(Rmax, Rn, Xn, R(j));
        else
            X(j) = 0.5;
        end
    end
    Y = V_HRev(Vm, Rmax, B, R, length(R), X);
end

function Y = velocityModel_p0(beta, X)

    [Y_Input, Vm, dirTrans, Vtr, dirShr, Vsh, sizeInput] = fetchModelInputs(X);
    [A_Trans, Phi_Trans] = getVtrParameters(beta, Vtr, sizeInput);
    
    if length(beta) > 4
        betaSh2 = beta(5:6); betaSh2_A = betaSh2(1); betaSh2_Phi = betaSh2(2);
        A_Shr2 = betaSh2_A(1)*ones(sizeInput); Phi_Shr2 = betaSh2_Phi(1)*ones(sizeInput);
    end
    
    if length(beta) > 6
        betaSh3 = beta(7:8); betaSh3_A = betaSh3(1); betaSh3_Phi = betaSh3(2);   
        A_Shr3 = betaSh3_A(1)*ones(sizeInput); Phi_Shr3 = betaSh3_Phi(1)*ones(sizeInput);
    end
    
    velTr = A_Trans .* cos(dirTrans - Phi_Trans);
    
    if length(beta) == 4
        velAsym = velTr;
    elseif length(beta) == 6
        velSh2 = A_Shr2 .* cos(2 * (dirShr - Phi_Shr2));
        velAsym = velTr + velSh2;
    elseif length(beta) == 8
        velSh2 = A_Shr2 .* cos(2 * (dirShr - Phi_Shr2));
        velSh3 = A_Shr3 .* cos(3 * (dirShr - Phi_Shr3));
        velAsym = velTr + velSh2 + velSh3;
    end
    
    Y = Y_Input .* (Vm + velAsym);
end

function Y = velocityModel_p1(beta, X)

    [Y_Input, Vm, dirTrans, Vtr, dirShr, Vsh, sizeInput] = fetchModelInputs(X);
    [A_Trans, Phi_Trans] = getVtrParameters(beta, Vtr, sizeInput);

    if length(beta) > 4
        betaSh1 = beta(5:8); betaSh1_A = betaSh1(1:2); betaSh1_Phi = betaSh1(3:4);
        A_Shr1 = betaSh1_A(1)*ones(sizeInput) + betaSh1_A(2)*Vsh; 
        Phi_Shr1 = betaSh1_Phi(1)*ones(sizeInput) + betaSh1_Phi(2)*Vsh;
        velSh1 = A_Shr1 .* cos(1 * (dirShr - Phi_Shr1));
    end    
    
    if length(beta) > 8
        betaSh2 = beta(9:12); betaSh2_A = betaSh2(1:2); betaSh2_Phi = betaSh2(3:4);
        A_Shr2 = betaSh2_A(1)*ones(sizeInput) + betaSh2_A(2)*Vsh; 
        Phi_Shr2 = betaSh2_Phi(1)*ones(sizeInput) + betaSh2_Phi(2)*Vsh;
        velSh2 = A_Shr2 .* cos(2 * (dirShr - Phi_Shr2));
    end
    
    if length(beta) > 12
        betaSh3 = beta(13:16); betaSh3_A = betaSh3(1:2); betaSh3_Phi = betaSh3(3:4);
        A_Shr3 = betaSh3_A(1)*ones(sizeInput) + betaSh3_A(2)*Vsh; 
        Phi_Shr3 = betaSh3_Phi(1)*ones(sizeInput) + betaSh3_Phi(2)*Vsh;
        velSh3 = A_Shr3 .* cos(3 * (dirShr - Phi_Shr3));
    end
    
    velTr = A_Trans .* cos(dirTrans - Phi_Trans); 
    if length(beta) == 4
        velAsym = velTr;
    elseif length(beta) == 8
        velAsym = velTr + velSh1;
    elseif length(beta) == 12
        velAsym = velTr + velSh1 + velSh2;
    elseif length(beta) == 16
        velAsym = velTr + velSh1 + velSh2 + velSh3;
    end
    
    Y = Y_Input .* (Vm + velAsym);
end

% function Y = velocityModel_p01(beta, X)
% 
%     [Y_Input, Vm, dirTrans, Vtr, dirShr, Vsh, sizeInput] = fetchModelInputs(X);
%     [A_Trans, Phi_Trans] = getVtrParameters(beta, Vtr, sizeInput);
%     
%     betaSh2 = beta(5:6); betaSh2_A = betaSh2(1); betaSh2_Phi = betaSh2(2);
%     A_Shr2 = betaSh2_A(1)*ones(sizeInput); Phi_Shr2 = betaSh2_Phi(1)*ones(sizeInput);
%     
%     betaSh3 = beta(7:10); betaSh3_A = betaSh3(1:2); betaSh3_Phi = betaSh3(3:4);      
%     A_Shr3 = betaSh3_A(1)*ones(sizeInput) + betaSh3_A(2)*Vsh; 
%     Phi_Shr3 = betaSh3_Phi(1)*ones(sizeInput) + betaSh3_Phi(2)*Vsh;
%     
%     velTr = A_Trans .* cos(dirTrans - Phi_Trans);
%     velSh2 = A_Shr2 .* cos(2 * (dirShr - Phi_Shr2));
%     velSh3 = A_Shr3 .* cos(3 * (dirShr - Phi_Shr3));
%     velAsym = velTr + velSh2 + velSh3;
%     
%     Y = Y_Input .* (Vm + velAsym);
% end
% 
% function Y = velocityModel_p10(beta, X)
% 
%     [Y_Input, Vm, dirTrans, Vtr, dirShr, Vsh, sizeInput] = fetchModelInputs(X);
%     [A_Trans, Phi_Trans] = getVtrParameters(beta, Vtr, sizeInput);
%     
%     betaSh2 = beta(5:8); betaSh2_A = betaSh2(1:2); betaSh2_Phi = betaSh2(3:4);
%     A_Shr2 = betaSh2_A(1)*ones(sizeInput) + betaSh2_A(2)*Vsh; 
%     Phi_Shr2 = betaSh2_Phi(1)*ones(sizeInput) + betaSh2_Phi(2)*Vsh;
%     
%     betaSh3 = beta(9:10); betaSh3_A = betaSh3(1); betaSh3_Phi = betaSh3(2);   
%     A_Shr3 = betaSh3_A(1)*ones(sizeInput); Phi_Shr3 = betaSh3_Phi(1)*ones(sizeInput);
%     
%     velTr = A_Trans .* cos(dirTrans - Phi_Trans);
%     velSh2 = A_Shr2 .* cos(2 * (dirShr - Phi_Shr2));
%     velSh3 = A_Shr3 .* cos(3 * (dirShr - Phi_Shr3));
%     velAsym = velTr + velSh2 + velSh3;
%     
%     Y = Y_Input .* (Vm + velAsym);
% end

function [Y_Input, Vm, dirTrans, Vtr, dirShr, Vsh, sizeInput] = fetchModelInputs(X)
    Y_Input = X(:,1); Vm = X(:,2); dirTrans = X(:,3); Vtr = X(:,4); 
    dirShr = X(:,5); Vsh = X(:,6);
    sizeInput = size(Y_Input);
end

function [A_Trans, Phi_Trans] = getVtrParameters(beta, Vtr, sizeInput)
    betaTr = beta(1:4); betaTr_A = betaTr(1:2); betaTr_Phi = betaTr(3:4);
    A_Trans = betaTr_A(1)*ones(sizeInput) + betaTr_A(2)*Vtr; 
    Phi_Trans = betaTr_Phi(1)*ones(sizeInput) + betaTr_Phi(2)*Vtr;
end

function models = getResidualStatisticsTraining(models, index, filePath)
    load(filePath);
    
    [models(index).MSE_Train, models(index).prctile16_Train, models(index).prctile50_Train, ...
        models(index).prctile84_Train, models(index).resRadialBin_Train] = ...
        getResidualStatistics(models(index).resTrain, StormTrain);

    [models(index).MSE_Val, models(index).prctile16_Val, models(index).prctile50_Val, ...
        models(index).prctile84_Val, models(index).resRadialBin_Val] = ...
        getResidualStatistics(models(index).resVal, StormVal);
end

function [MSE, prctile16, prctile50, prctile84, resRadialBin] = ...
        getResidualStatistics(residual, Storm)
    MSE = residual.^2/length(residual);
    prctile16 = prctile(residual.^2, 16);
    prctile50 = prctile(residual.^2, 50);
    prctile84 = prctile(residual.^2, 84);
    
    RNorm = getNormalizedRadius(Storm);

    idx1 = find(RNorm < 1);
    idx2 = find(RNorm >= 1 & RNorm < 2);
    idx3 = find(RNorm >= 2 & RNorm < 4);
    idx4 = find(RNorm >= 4 & RNorm < 7);  

    sqResidual = residual.^2;

    resRadialBin = [mean(sqResidual(idx1)); mean(sqResidual(idx2)); ...
        mean(sqResidual(idx3)); mean(sqResidual(idx4))]; 
    resRadialBin = [0; mean(sqResidual(idx2)); ...
        mean(sqResidual(idx3)); mean(sqResidual(idx4))];    
end

function RNorm = getNormalizedRadius(Storm)
    Rmax = []; R = [];
    [Rcut, ~, ~, ~, ~, ~] = CutoffCalc2(Storm);
    for i = 1:length(Storm)
        idx = find(Storm(i).AsymR < Rcut(i));
        Rmax = [Rmax; Storm(i).Rmax*ones(length(idx), 1)];
        R = [R; Storm(i).AsymR(idx)'];
    end
    RNorm = R./Rmax;
end