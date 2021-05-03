function models = ResultsTrainTestVecAdd(x0, b_l, b_u, noWvn, opt, Rbound)
    
    tic
    models = getResultsTrain(x0, b_l, b_u, noWvn, opt, Rbound);
    disp('Get Results')
    
    models = ResultsTestHWind(models, noWvn, opt, Rbound);
    disp('Testing HWind')
    
%     [MSE_KatrinaWRF, MSE_SandyWRF] = ResultsTestWRF(models(6).x, noWvn, opt);
%     models(6).MSE_KatrinaWRF = MSE_KatrinaWRF;
%     models(6).MSE_SandyWRF = MSE_SandyWRF;
%     disp('Testing WRF')
    toc
end

function models = getResultsTrain(x0, b_l, b_u, noWvn, opt, Rbound)
    idx = 1:1:5;
    path = '/Users/changd/Dropbox (MIT)/2019_JAMC/MATLAB_revision/';
    models = struct([]);
    
    for i = 1:5
%         i
        filePath = [path, 'Training Data/TrainingSets/Training', ...
            num2str(i), '_', opt, '.mat'];
        [models(i).x, ci, models(i).res, models(i).resVal, models(i).AIC, ...
            models(i).MSE_train, models(i).MSE_val] = ...
            getResults(filePath, x0(:,i), b_l, b_u, noWvn, Rbound);
        models(i).ci = (ci(:,2) - ci(:,1))/2;
        disp(['Done with training set ', num2str(i)])
%         disp(['Parameters: ', models(i).x])
    end
    
    filePath = [path, 'Training Data/TrainingSets/TrainingAll', opt, '.mat'];
    [models(6).x, ci, models(6).res, models(6).AIC, models(6).MSE_train] = ...
    getResultsTrainingAll(filePath, x0(:,6), b_l, b_u, noWvn, Rbound);
    models(6).ci = (ci(:,2) - ci(:,1))/2;
    disp(['Done with full training set '])
end

function [x, ci, rTrain, AIC, MSE_train] = getResultsTrainingAll ...
        (filePath, x0, b_l, b_u, noWvn, Rbound)

    load(filePath); lambda=0;
    
    % Get inputs for optimization problem (training set)    
    [X, V, R] = getInputs(Storm, Rbound); 
    V_Hol = updateStormAsymTrans(Storm, Rbound); N_Train = length(V);
    V = V-V_Hol;
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [A, b] = getConstraints(Storm);

    % Solve optimization problem
    velocityModel(x0); velocityModelHessian(x0,lambda);
    options = optimoptions('fmincon','Algorithm','interior-point','Hessian','user-supplied',...
        'GradObj','on','Display','iter','HessFcn',@velocityModelHessian,'MaxIter',10000);
    problem = createOptimProblem('fmincon',...
        'lb',b_l,'ub',b_u,'x0',x0,'Aineq',A,'bineq',b,'options',options);
    [x, fTrain, ~, ~, ~ , ~ , ~] = fmincon(@velocityModel, x0, A, b, [], [], b_l, ...
        b_u, [], options);
    
    % Get confidence interval / AIC
    [X, V, R] = getInputs(Storm, Rbound); N_Train = length(V);
    V_Hol = updateStormAsymTrans(Storm, Rbound); 
    V = V-V_Hol;
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    velocityModel(x);
    load('residuals.mat'); rTrain = r;
    
    if sum(abs(x(5:end))) > 0
        ci = nlparci(x, r, 'jacobian', J); k = length(x) + 4;
    else
        ci = nlparci(x(1:4), r, 'jacobian', J(:,1:4)); k = 8;
    end
    AIC = N_Train*log(2*fTrain/N_Train) + 2*k;
    
    % Get SSE
    MSE_train = 2*fTrain/N_Train; 
end

function [x, ci, rTrain, rVal, AIC, MSE_train, MSE_val] = getResults(filePath, x0, b_l, b_u, noWvn, Rbound)

    load(filePath); lambda=0;
    
    % Get inputs for optimization problem (training set)    
    [X, V, R] = getInputs(StormTrain, Rbound); 
    V_Hol = updateStormAsymTrans(StormTrain, Rbound); N_Train = length(V);
    V = V-V_Hol;
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [A, b] = getConstraints(StormTrain);

    % Solve optimization problem
    velocityModel(x0); velocityModelHessian(x0,lambda);
    options = optimoptions('fmincon','Algorithm','interior-point','Hessian','user-supplied',...
        'GradObj','on','Display','iter','HessFcn',@velocityModelHessian,'MaxIter',10000);
    problem = createOptimProblem('fmincon',...
        'lb',b_l,'ub',b_u,'x0',x0,'Aineq',A,'bineq',b,'options',options);
    [x, fTrain, ~, ~, ~ , ~ , ~] = fmincon(@velocityModel, x0, A, b, [], [], b_l, ...
        b_u, [], options);
    
    % Get confidence interval / AIC
    [X, V, R] = getInputs(StormTrain, Rbound); N_Train = length(V);
    V_Hol = updateStormAsymTrans(StormTrain, Rbound); 
    V = V-V_Hol;
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    velocityModel(x);
    load('residuals.mat'); rTrain = r;
    
    SSE = 2*fTrain/N_Train; 
    if sum(abs(x(5:end))) > 0
        ci = nlparci(x, r, 'jacobian', J); k = length(x) + 4;
    else
        ci = nlparci(x(1:4), r, 'jacobian', J(:,1:4)); k = 8;
    end
    AIC = N_Train*log(2*fTrain/N_Train) + 2*k;
%     AIC = 2*fTrain/var(r) + 2*k;
    
    % Get results (validation)
    [X, V, R] = getInputs(StormVal, Rbound); 
    V_Hol = updateStormAsymTrans(StormVal, Rbound); N_Val = length(V);
    V = V-V_Hol;
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [fVal, ~] = velocityModel(x); 
    load('residuals.mat'); rVal = r;
    
    % Get SSE
    MSE_train = 2*fTrain/N_Train; MSE_val = 2*fVal/N_Val;
end

function models = ResultsTestHWind(models, noWvn, opt, RBound)

    x = zeros(length(models(1).x), length(models));
    for i = 1:length(models)
        x(:,i) = models(i).x;
    end

    path = '/Users/changd/Dropbox (MIT)/2019_JAMC/MATLAB_revision/';
    
    i = 6;
    load([path 'Training Data/TestStorms/IreneHW_' opt '.mat']);
    models(i).SSE_Irene = ResultsTest(Storm, x(:,i), noWvn, RBound);
    for j = 1:length(Storm)
        models(i).SSE_IreneSnapshot(j) = ResultsTest(Storm(j), x(:,i), noWvn, RBound);
        load('asymParamTr.mat'); load('asymParamSh.mat');
        models(i).ATrans_Irene(j) = unique(A_Trans);
        models(i).PhiTrans_Irene(j) = unique(Phi_Trans);
        models(i).AShr_Irene(j) = unique(A_Shr1);
        models(i).PhiShr_Irene(j) = unique(Phi_Shr1);
    end
    
    load([path 'Training Data/TestStorms/IsabelHW_SHIPS.mat']);
    models(i).SSE_Isabel = ResultsTest(Storm, x(:,i), noWvn, RBound);
    for j = 1:length(Storm)
        models(i).SSE_IsabelSnapshot(j) = ResultsTest(Storm(j), x(:,i), noWvn, RBound);
        load('asymParamTr.mat'); load('asymParamSh.mat');
        models(i).ATrans_Isabel(j) = unique(A_Trans);
        models(i).PhiTrans_Isabel(j) = unique(Phi_Trans);
        models(i).AShr_Isabel(j) = unique(A_Shr1);
        models(i).PhiShr_Isabel(j) = unique(Phi_Shr1);
    end        
    
    load([path 'Training Data/TestStorms/KatrinaHW_' opt '.mat']);
    models(i).SSE_Katrina = ResultsTest(Storm, x(:,i), noWvn, RBound);
    for j = 1:length(Storm)
        models(i).SSE_KatrinaSnapshot(j) = ResultsTest(Storm(j), x(:,i), noWvn, RBound);
        load('asymParamTr.mat'); load('asymParamSh.mat');
        models(i).ATrans_Katrina(j) = unique(A_Trans);
        models(i).PhiTrans_Katrina(j) = unique(Phi_Trans);
        models(i).AShr_Katrina(j) = unique(A_Shr1);
        models(i).PhiShr_Katrina(j) = unique(Phi_Shr1);
    end
    
    load([path 'Training Data/TestStorms/SandyHW_' opt '.mat']);
    models(i).SSE_Sandy = ResultsTest(Storm, x(:,i), noWvn, RBound);
    load([path 'Training Data/TestStorms/SandyHWTrain' '.mat']);
    for j = 1:length(Storm)
        models(i).SSE_SandySnapshot(j) = ResultsTest(Storm(j), x(:,i), noWvn, RBound);
        load('asymParamTr.mat'); load('asymParamSh.mat');
        models(i).ATrans_Sandy(j) = unique(A_Trans);
        models(i).PhiTrans_Sandy(j) = unique(Phi_Trans);
        models(i).AShr_Sandy(j) = unique(A_Shr1);
        models(i).PhiShr_Sandy(j) = unique(Phi_Shr1);
    end
end

function [SSE_Katrina, SSE_Sandy] = ResultsTestWRF(x, noWvn, opt)

    path = '/Users/changd/Dropbox (MIT)/2019_JAMC/MATLAB_revision/';
    load([path 'Training Data/TestStorms/KatrinaWRF_' opt '.mat']);
    SSE_Sandy = zeros(10,1);
    idxSandy = [13 14 16 30 34 42 46 52 56 59];
    
    SSE_Katrina = ResultsTestWRF_Iter(Storm, x, noWvn, [0 Inf]);
    for j = 1:10
%         j
        load([path 'Training Data/TestStorms/TestSandyWRF_' opt '_' num2str(idxSandy(j)) '.mat']);
        SSE_Sandy(j) = ResultsTestWRF_Iter(Storm, x, noWvn, [0 Inf]);
    end
end

function [SSE] = ResultsTest(Storm, xOpt, noWvn, Rbound)
    [X, V, R] = getInputs(Storm, Rbound); 
    V_Hol = updateStormAsymTrans(Storm, Rbound);
    V = V-V_Hol;
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [f, ~] = velocityModel(xOpt); 
    SSE = 2*f/length(V);
end

function [MSE, r] = ResultsTestWRF_Iter(Storm, xOpt, noWvn, Rbound)
    [X, V, R] = getInputs(Storm, Rbound);
    for i = 1:length(V)
        if V(i) > 23
            gf = 1 + 2.5*10e-10*(V(i)-23)^4.7;
            V(i) = V(i)/gf;
        end
    end
    
    V_Hol = updateStormAsymTrans(Storm, Rbound);
    V = V-V_Hol;
 
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [f, ~] = velocityModel(xOpt); 
    MSE = 2*f/length(V);
    load('residuals.mat');
end

function [X, V, R] = getInputs(Storm, Rbound)
    R = [];
    Y = []; Vm = []; dirTrans = []; Vtr = []; dirShr = []; Vsh = []; V = [];
    [Rcut, ~, ~, ~, ~, ~] = CutoffCalc2(Storm);
    for i = 1:length(Storm)
        [YIter, ~, idx] = VHolCNLS_Iter(Storm(i).AsymR, Storm(i).Vm, ...
            Storm(i).Rmax, Storm(i).B, Storm(i).Rn, Storm(i).Xn, Rcut(i), ...
            Rbound);

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

function [Y, R, idx] = VHolCNLS_Iter(R, Vm, Rmax, B, Rn, Xn, Rcut, Rbound)
    
    Xcalc = @(Rm, Rn, Xn, R) 0.5+(R-Rm)*(Xn-0.5)/(Rn-Rm);
    V_HRev=@(Vm,Rm,B,R,N,X) ((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^X;
    
    idx=find(R < Rcut & R/Rmax >= Rbound(1) & R/Rmax <= Rbound(2)); 
    R = R(idx);
%     disp(['Number of data points: ', num2str(length(R))])

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

function V = updateStormAsymTrans(Storm, Rbound)
    
    residuals = []; N = 0;
    [Rcut, ~, ~, ~, ~, ~] = CutoffCalc2(Storm);
    
    V = [];
    for i = 1:length(Storm)
        % Get lambdas with respect to north + Holland velocities 
        dir = getDirNormal(Storm, i);
        [Y, ~, idx] = VHolCNLS_Iter(Storm(i).AsymR, Storm(i).Vm, Storm(i).Rmax, ...
            Storm(i).B, Storm(i).Rn, Storm(i).Xn, Rcut(i), Rbound);
        dir = dir(idx);
        
        % Break down MF to u- and v-components
        [u, v] = getMF_Components(Y*Storm(i).Vm, dir);  
        
        % Break down translation vector to u- and v-components
        [uTrans, vTrans] = getTranslationComponents(Storm(i).Vtrans, ...
            Storm(i).Vtransdir);
    
        % Vector addition
        uAsym = u + uTrans*ones(size(u));
        vAsym = v + vTrans*ones(size(u));
        
        % Get velocity magnitude
        velAsym = sqrt(uAsym.^2 + vAsym.^2);
        V = [V; velAsym];
    end
end

function [f, df] = velocityModel(beta)

    load('Input.mat');
    V_data = V;
    [Y_Input, Vm, dirTrans, Vtr, dirShr, Vsh, sizeInput] = fetchModelInputs(X);
    [A_Trans, Phi_Trans] = getVtrParameters(beta, Vtr, sizeInput);
    save('asymParamTr.mat', 'A_Trans', 'Phi_Trans');

    if length(beta) > 4
        betaSh1 = beta(5:8); betaSh1_A = betaSh1(1:2); betaSh1_Phi = betaSh1(3:4);
        A_Shr1 = betaSh1_A(1)*ones(sizeInput) + betaSh1_A(2)*Vsh; 
        Phi_Shr1 = betaSh1_Phi(1)*ones(sizeInput) + betaSh1_Phi(2)*Vsh;
        velSh1 = A_Shr1 .* cos(1 * (dirShr - Phi_Shr1));
        save('asymParamSh.mat', 'A_Shr1', 'Phi_Shr1');
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
    
%     V_pred = Y_Input .* (Vm + velAsym);
    V_pred = Y_Input .* (velAsym);
    
    r = V_data - V_pred;
    f = 0.5*(r'*r);
    
    N = length(r);
    J_tr_A = zeros(N,2); J_tr_Phi = zeros(N,2);
    for i = 1:N
        J_tr_A(i,:) = -Y_Input(i) * cos(dirTrans(i) - Phi_Trans(i)) * [1 Vtr(i)];
        J_tr_Phi(i,:) = -Y_Input(i) * A_Trans(i) * sin(dirTrans(i) - Phi_Trans(i)) ...
            * [1 Vtr(i)];
    end
    J = [J_tr_A J_tr_Phi];
    
    if length(beta) > 4
        J_sh1_A = zeros(N,2); J_sh1_Phi = zeros(N,2);
        for i = 1:N
            J_sh1_A(i,:) = -Y_Input(i) * cos(1*(dirShr(i) - Phi_Shr1(i))) * [1 Vsh(i)];
            J_sh1_Phi(i,:) = -Y_Input(i) * A_Shr1(i) * sin(1*(dirShr(i) - Phi_Shr1(i))) ...
            * [1 Vsh(i)];
        end
        J = [J J_sh1_A J_sh1_Phi];
        save('InputHessian.mat', 'r', 'J', 'Y_Input', 'A_Trans', 'Phi_Trans', ...
            'A_Shr1', 'Phi_Shr1', 'dirTrans', 'dirShr', 'Vtr', 'Vsh');
    end
    
    if length(beta) > 8
        J_sh2_A = zeros(N,2); J_sh2_Phi = zeros(N,2);
        for i = 1:N
            J_sh2_A(i,:) = -Y_Input(i) * cos(2*(dirShr(i) - Phi_Shr2(i))) * [1 Vsh(i)];
            J_sh2_Phi(i,:) = -Y_Input(i) * 2 * A_Shr2(i) * sin(2*(dirShr(i) - Phi_Shr2(i))) ...
            * [1 Vsh(i)];
        end
        J = [J J_sh2_A J_sh2_Phi];
        save('InputHessian.mat', 'r', 'J', 'Y_Input', 'A_Trans', 'Phi_Trans', ...
            'A_Shr1', 'Phi_Shr1', 'A_Shr2', 'Phi_Shr2', 'dirTrans', 'dirShr', ...
            'Vtr', 'Vsh');
    end
    
    if length(beta) > 12
        J_sh3_A = zeros(N,2); J_sh3_Phi = zeros(N,2);
        for i = 1:N
            J_sh3_A(i,:) = -Y_Input(i) * cos(3*(dirShr(i) - Phi_Shr3(i))) * [1 Vsh(i)];
            J_sh3_Phi(i,:) = -Y_Input(i) * 3 * A_Shr3(i) * sin(3*(dirShr(i) - Phi_Shr3(i))) ...
            * [1 Vsh(i)];
        end
        J = [J J_sh3_A J_sh3_Phi];
        save('InputHessian.mat', 'r', 'J', 'Y_Input', 'A_Trans', 'Phi_Trans', ...
            'A_Shr1', 'Phi_Shr1', 'A_Shr2', 'Phi_Shr2', 'A_Shr3', 'Phi_Shr3', ...
            'dirTrans', 'dirShr', 'Vtr', 'Vsh');        
    end
    
    df = J'*r;
    
    save('residuals.mat', 'r', 'J');
end

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

function dir = getDirNormal(Storm, index)
    dir = Storm(index).AsymDir*(360/(2*pi));
    dir = dir + Storm(index).Vtransdir*ones(size(Storm(index).AsymDir));
    idx = find(dir < 0);
    dir(idx) = dir(idx) + 360*ones(size(dir(idx)));
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