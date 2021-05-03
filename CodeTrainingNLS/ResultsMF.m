function [residuals, MSE_Hwind, MSE_WRF] = ResultsMF(opt)
    residuals = getResultsTrain(opt);
    
    [MSE_Irene, MSE_KatrinaHW, MSE_SandyHW, MSE_Isabel] = ResultsTestHWind(opt);
    MSE_Hwind = [MSE_Irene, MSE_KatrinaHW, MSE_SandyHW, MSE_Isabel];
    
%     [MSE_KatrinaWRF, MSE_SandyWRF] = ResultsTestWRF(opt);
    MSE_KatrinaWRF = 0; MSE_SandyWRF = 0;
    MSE_WRF = [MSE_KatrinaWRF, MSE_SandyWRF];
end

function residuals = getResultsTrain(opt)
    path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';

%     idx = 1:5;
%     MSE_train = zeros(1, length(idx)); MSE_val = zeros(1, length(idx));
%     for i = 1:length(idx)
%         i
        filePath = [path, 'Training Data/TrainingSets/TrainingAllGFS.mat'];
        load(filePath);
        [~, ~, residuals] = updateStormAsymTrans(Storm);
%         [MSE_val(i), SSE_val, residuals_val] = updateStormAsymTrans(StormVal);
%     end
        disp('check')
end

function [MSE_Irene, MSE_Katrina, MSE_Sandy, MSE_Isabel] ...
        = ResultsTestHWind(opt)

    MSE_Irene = zeros(1,5);
    MSE_Katrina = zeros(1,5);
    MSE_Sandy = zeros(1,5);
    MSE_Isabel = zeros(1,5);

    path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
    
    for i = 1:5
        load([path 'Training Data/TestStorms/IreneHW_' opt '.mat']);
        [MSE_Irene(i), ~, ~] = updateStormAsymTrans(Storm);
        load([path 'Training Data/TestStorms/KatrinaHW_' opt '.mat']);
        [MSE_Katrina(i), ~, ~] = updateStormAsymTrans(Storm);
        load([path 'Training Data/TestStorms/SandyHW_' opt '.mat']);
        [MSE_Sandy(i), ~, ~] = updateStormAsymTrans(Storm);
        load([path 'Training Data/TestStorms/IsabelHW_SHIPS.mat']);
        [MSE_Isabel(i), ~, ~] = updateStormAsymTrans(Storm);
    end
end

function [MSE_Katrina, MSE_Sandy] = ResultsTestWRF(opt)

    MSE_Katrina = zeros(5,1);
    MSE_Sandy = zeros(5,10);
    idxSandy = [13 14 16 30 34 42 46 52 56 59];
    for i = 1:5
        %i
        load(['Training Data/KatrinaWRF_' opt '.mat']);
        [MSE_Katrina(i), ~, ~] = updateStormAsymTrans(Storm);
        for j = 1:10
            %j
            load(['Training Data/TestSandyWRF_' opt '_' num2str(idxSandy(j)) '.mat']);
            [MSE_Sandy(i,j), ~, ~] = updateStormAsymTrans(Storm);
        end
    end
end

function [MSE, SSE, residuals] = updateStormAsymTrans(Storm)
    
    residuals = []; N = 0;
    [Rcut, ~, ~, ~, ~, ~] = CutoffCalc2(Storm);
    for i = 1:length(Storm)
        % Get lambdas with respect to north + Holland velocities 
        dir = getDirNormal(Storm, i);
        [Y, ~, idx] = VHolCNLS_Iter(Storm(i).AsymR, Storm(i).Vm, Storm(i).Rmax, ...
            Storm(i).B, Storm(i).Rn, Storm(i).Xn, Rcut(i));
        dir = dir(idx);
        
        residuals_Iter = Y' - Storm(i).AsymV(idx)';
        residuals = [residuals; residuals_Iter]; N = N + length(idx);
    end
    
    SSE = sum(residuals);
    MSE = SSE/N;
end

function dir = getDirNormal(Storm, index)
    dir = Storm(index).AsymDir*(360/(2*pi));
    dir = dir + Storm(index).Vtransdir*ones(size(Storm(index).AsymDir));
    idx = find(dir < 0);
    dir(idx) = dir(idx) + 360*ones(size(dir(idx)));
end

function [Y, R, idx] = VHolCNLS_Iter(R, Vm, Rmax, B, Rn, Xn, Rcut)
    
    Xcalc = @(Rm, Rn, Xn, R) 0.5+(R-Rm)*(Xn-0.5)/(Rn-Rm);
    V_HRev=@(Vm,Rm,B,R,N,X) Vm*((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^X;
    
    idx=find(R < Rcut); R = R(idx);

    X = zeros(1, length(R));
    for j = 1:length(R)
        if R(j) > Rmax;
            X(j) = Xcalc(Rmax, Rn, Xn, R(j));
        else
            X(j) = 0.5;
        end
    end
    Y = V_HRev(Vm, Rmax, B, R, length(R), X);
end

function [u, v] = getMF_Components(V, dir)
    
    u = zeros(length(V), 1); v = zeros(length(V), 1);
    for j = 1:length(dir)
        dirRotation = dir(j) + 270;
        if dirRotation > 360
            dirRotation = dirRotation - 360;
        end
        if dir(j) < 90
            angle = dirRotation;
            u(j) = V(j)*sind(angle);
            v(j) = V(j)*cosd(angle);
        elseif dir(j) >= 90 && dir(j) < 180
            angle = 180 - dirRotation;
            u(j) = V(j)*sind(angle);
            v(j) = -V(j)*cosd(angle);
        elseif dir(j) >= 180 && dir(j) < 270
            angle = dirRotation - 180;
            u(j) = -V(j)*sind(angle);
            v(j) = -V(j)*cosd(angle);
        elseif dir(j) >= 270 && dir(j) <= 360
            angle = 360 - dirRotation;
            u(j) = -V(j)*sind(angle);
            v(j) = V(j)*cosd(angle);
        else
            disp('error')
        end
    end
end

function [uTrans, vTrans] = getTranslationComponents(Vtrans, Vtransdir)

    if Vtransdir < 90
        angle = Vtransdir;
        uTrans = Vtrans*sind(angle);
        vTrans = Vtrans*cosd(angle);
    elseif Vtransdir >= 90 && Vtransdir < 180
        angle = 180 - Vtransdir;
        uTrans = Vtrans*sind(angle);
        vTrans = -Vtrans*cosd(angle);
    elseif Vtransdir >= 180 && Vtransdir < 270
        angle = Vtransdir - 180;
        uTrans = -Vtrans*sind(angle);
        vTrans = -Vtrans*cosd(angle);
    elseif Vtransdir >= 270 && Vtransdir <= 360
        angle = 360 - Vtransdir;
        uTrans = -Vtrans*sind(angle);
        vTrans = Vtrans*cosd(angle);
    end
end