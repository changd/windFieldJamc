function models = ResultsTrainTest(x0, b_l, b_u, noWvn, opt, Rbound)
    
    tic
    models = getResultsTrain(x0, b_l, b_u, noWvn, opt, Rbound);
    disp('Get Results')

    
    models = ResultsTestHWind(models, noWvn, opt, Rbound);
    disp('Testing HWind')
    
    [MSE_KatrinaWRF, MSE_SandyWRF, meanResKatrinaWRF, stdResKatrinaWRF, ...
        meanResSandyWRF, stdResSandyWRF] = ResultsTestWRF(models(6).x, noWvn, opt);
    disp('Testing WRF')
    models(6).MSE_KatrinaWRF = MSE_KatrinaWRF;
    models(6).MSE_SandyWRF = MSE_SandyWRF;
    models(6).meanResKatrinaWRF = meanResKatrinaWRF;
    models(6).stdResKatrinaWRF = stdResKatrinaWRF;
    models(6).meanResSandyWRF = meanResSandyWRF;
    models(6).stdResSandyWRF = stdResSandyWRF;
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
    [X, V, R] = getInputs(Storm, Rbound); N_Train = length(V);
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
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    velocityModel(x);
    load('residuals.mat'); rTrain = r;
    
    if sum(abs(x(5:end))) > 0
        ci = nlparci(x, r, 'jacobian', J); k = length(x) + 4;
    else
        ci = nlparci(x(1:4), r, 'jacobian', J(:,1:4)); k = 8;
    end
    AIC = N_Train*log(2*fTrain/N_Train) + 2*k;
    
    % Get MSE
    MSE_train = 2*fTrain/N_Train;
end

function [x, ci, rTrain, rVal, AIC, MSE_train, MSE_val] = getResults(filePath, x0, b_l, b_u, noWvn, Rbound)

    load(filePath); lambda=0; 
    
    % Get inputs for optimization problem (training set)    
    [X, V, R] = getInputs(StormTrain, Rbound); N_Train = length(V);
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
    [X, V, R] = getInputs(StormVal, Rbound); N_Val = length(V);
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [fVal, ~] = velocityModel(x); 
    load('residuals.mat'); rVal = r;
    
    % Get MSE
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
end

function [SSE_Katrina, SSE_Sandy, meanResKatrina, stdResKatrina, ...
        meanResSandy, stdResSandy] = ResultsTestWRF(x, noWvn, opt)

    path = '/Users/changd/Dropbox (MIT)/2019_JAMC/MATLAB_revision/';
    SSE_Sandy = zeros(10,1);
    idxSandy = [13 14 16 30 34 42 46 52 56 59];
    %i
    load([path 'Training Data/TestStorms/KatrinaWRF_' opt '.mat']);
    [SSE_Katrina, res] = ResultsTestWRF_Iter(Storm, x, noWvn, [0 Inf]);
    meanResKatrina = mean(res); stdResKatrina = std(res);
    meanResSandy = zeros(10,1); stdResSandy = zeros(10,1);
    for j = 1:10
%         j
        load([path 'Training Data/TestStorms/TestSandyWRF_' opt '_' num2str(idxSandy(j)) '.mat']);
        [SSE_Sandy(j), res] = ResultsTestWRF_Iter(Storm, x, noWvn, [0 Inf]);
        meanResSandy(j) = mean(res); stdResSandy(j) = std(res);
    end
end

function [MSE, r] = ResultsTest(Storm, xOpt, noWvn, Rbound)
    [X, V, R] = getInputs(Storm, Rbound); 
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [f, ~] = velocityModel(xOpt); 
    MSE = 2*f/length(V);
    load('residuals.mat');
end

function [MSE, r] = ResultsTestWRF_Iter(Storm, xOpt, noWvn, Rbound)
    [~,L]=size(xOpt);
    SSE=zeros(1,L);

    [X, V, R] = getInputs(Storm, Rbound);

    for j = 1:length(V)
        if V(j) > 23
            gf = 1 + 2.5*10e-10*(V(j)-23)^4.7;
            V(j) = V(j)/gf;
        end
    end

    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [f, ~] = velocityModel(xOpt); 
    MSE = 2*f/length(V);
    load('residuals.mat');
end

function [X, V, R] = getInputs(Storm, Rbound)
    R = [];
    Y = []; Vm = []; dirTrans = []; Vtr = []; dirShr = []; Vsh = []; V = [];
    [Rcut, ~, ~, ~, ~, ~] = CutoffCalc2(Storm);
%     Storm = LocalAverage_Rev(Storm);
   
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
%         V_iter = vectorSubtraction(Storm, i);
        V = [V; Storm(i).AsymV(idx)'];
        R = [R; Storm(i).AsymR(idx)'/Storm(i).Rmax];
    end
    X = [Y Vm dirTrans Vtr dirShr Vsh];
end

function V_iter = vectorSubtraction(Storm, index, dataIndices)
    V = Storm(index).AsymV;
    dir = Storm(index).AsymTheta;
    
    % Break down MF to u- and v-components
    [u, v] = getVelocityDataComponents(V, dir);     
    
    % Break down translation vector to u- and v-components
    [uTrans, vTrans] = getTranslationComponents(Storm(index).Vtrans, ...
        Storm(index).Vtransdir);    
    
    % Vector addition
    uAsym = u + uTrans*ones(size(u));
    vAsym = v + vTrans*ones(size(u));

    % Get velocity magnitude
    velAsym = sqrt(uAsym.^2 + vAsym.^2);
    
    V_iter = velAsym;
end

function [u, v] = getVelocityDataComponents(V, dir)
    u = zeros(length(V), 1); v = zeros(length(V), 1);
    for i = 1:length(V)
        if dir(i) < 90
            angle = dir(i);
            u(i) = -V(i).*cos(angle);
            v(i) = V(i).*sin(angle);
        elseif dir(i) >= 90 && dir(i) < 180
            angle = 180 - dir(i);
            u(i) = V(i).*cos(angle);
            v(i) = V(i).*sin(angle);
        elseif dir(i) >= 180 && dir(i) < 270
            angle = dir(i) - 180;
            u(i) = V(i).*cos(angle);
            v(i) = -V(i).*sin(angle);
        elseif dir(i) >= 270 && dir(i) <= 360
            angle = 360 - dir(i);
            u(i) = -V(i).*cos(angle);
            v(i) = -V(i).*sin(angle);
        end
    end
end

function [uTrans, vTrans] = getTranslationComponents(Vtrans, Vtransdir)

    if Vtransdir < 90
        angle = Vtransdir;
        uTrans = Vtrans.*sind(angle);
        vTrans = Vtrans.*cosd(angle);
    elseif Vtransdir >= 90 && Vtransdir < 180
        angle = 180 - Vtransdir;
        uTrans = Vtrans.*sind(angle);
        vTrans = -Vtrans.*cosd(angle);
    elseif Vtransdir >= 180 && Vtransdir < 270
        angle = Vtransdir - 180;
        uTrans = -Vtrans.*sind(angle);
        vTrans = -Vtrans.*cosd(angle);
    elseif Vtransdir >= 270 && Vtransdir <= 360
        angle = 360 - Vtransdir;
        uTrans = -Vtrans.*sind(angle);
        vTrans = Vtrans.*cosd(angle);
    end
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

function Storm = LocalAverage_Rev(Storm)
    % Calculate locally averaged R, Rst, Dir, V, Vavg, Ravg

    divDir = 33;
    lengths=zeros(length(Storm),1);
    for i=1:length(Storm)
    %     i
        R=[0:5:35 40:(600-40)/25:600];
        dir=-180:360/divDir:180;

        V=zeros(length(R)-1,length(dir)-1); Rs=zeros(length(R)-1,length(dir)-1);
        dirs=zeros(length(R)-1,length(dir)-1); dirShrs=zeros(length(R)-1,length(dir)-1);
        thetas = zeros(length(R)-1, length(dir)-1); index=zeros(length(R)-1,length(dir)-1);
        for j=1:length(R)-1
            for k=1:length(dir)-1
                idx=find(Storm(i).R>=R(j) & Storm(i).R<R(j+1) & Storm(i).dirTrans>=dir(k) & Storm(i).dirTrans<dir(k+1));
                if isempty(idx)==0
                    V(j,k)=mean(Storm(i).V(idx));
                    Rs(j,k)=mean(Storm(i).R(idx));
                    dirs(j,k)=mean(Storm(i).dirTrans(idx));
                    dirShrs(j,k)=mean(Storm(i).dirShr(idx));
                    thetas(j,k) = mean(Storm(i).dir(idx));
                end
            end
            idx=find(V(j,:)~=0);
            if length(idx)==0
                Vavg(j)=0;
            else          
                Vavg(j)=mean(V(j,idx));
            end
            index(j,:)=j*ones(1,length(dir)-1);
        end

        [lengthR,lengthV]=size(V);
        Vrep=reshape(V,[lengthR*lengthV,1]); Rrep=reshape(Rs,[lengthR*lengthV,1]);
        dirrep=reshape(dirs,[lengthR*lengthV,1]); dirShrrep=reshape(dirShrs,[lengthR*lengthV,1]);
        thetasrep = reshape(thetas, [lengthR*lengthV, 1]);
        
        idx=find(Vrep~=0); lengths(i)=length(idx);
        Storm(i).AsymRst=Rrep(idx)'/Storm(i).Rmax;
        Storm(i).AsymR=Rrep(idx)'+0.1*ones(size(idx'));
        Storm(i).AsymDir=dirrep(idx)'*(pi/180);
        Storm(i).AsymDirShr=dirShrrep(idx)'*(pi/180);
        Storm(i).AsymV=Vrep(idx)';
        Storm(i).AsymTheta = thetasrep(idx)'*(pi/180);
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
    
    V_pred = Y_Input .* (Vm + velAsym);
    
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