function [x, SSE_train, SSE_val, SSE_IreneHW, SSE_KatrinaHW, SSE_SandyHW, ...
        SSE_IsabelHW, SSE_KatrinaWRF, SSE_SandyWRF] = ResultsTrainTest(x0, b_l, b_u, noWvn, opt)
    
    [x, SSE_train, SSE_val] = getResultsTrain(x0, b_l, b_u, noWvn, opt);
    disp('Get Results')

    %
    [SSE_IreneHW, SSE_KatrinaHW, SSE_SandyHW, SSE_IsabelHW] ...
        = ResultsTestHWind(x, noWvn, opt);
    disp('Testing HWind')
    
%     [SSE_KatrinaWRF, SSE_SandyWRF] = ResultsTestWRF(x, noWvn, opt);
%     disp('Testing WRF')
    SSE_KatrinaWRF = 0; SSE_SandyWRF = 0;
end

function [x, SSE_train, SSE_val] = getResultsTrain(x0, b_l, b_u, noWvn, opt)
    idx = 1:1:5;
    x = zeros(4+4*noWvn,length(idx));
    SSE_train = zeros(1,length(idx));
    SSE_val = zeros(1,length(idx));
    path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
    
    for i = 1:5
%         i
        filePath = [path, 'Training Data/TrainingSets/Training', ...
            num2str(i), '_', opt, '.mat'];
        [x(:,i), SSE_train(i), SSE_val(i)] = getResults(filePath, x0, b_l, b_u, noWvn);
    end
end

function [x, SSE_train, SSE_val] = getResults(filePath, x0, b_l, b_u, noWvn)

    load(filePath); lambda=0;
    
    % Get inputs for optimization problem (training set)
    [Rcut,~,~,~,~,~] = CutoffCalc2(StormTrain);
    [Y, ~] = VHolCNLS(StormTrain, Rcut);
    [D, DAsymC, Vm, V, StormTrain] = OptDataAsym1(StormTrain, Rcut);
    [A, b] = OptDataAsym2(DAsymC, noWvn);
    save('Input.mat', 'Y', 'D', 'Vm', 'V', 'noWvn'); [N_Train,~] = size(D);

    % Solve optimization problem
    ObjFnIPAsym(x0); ObjFnHessianIPAsym(x0,lambda);
    options = optimoptions('fmincon','Algorithm','interior-point','Hessian','user-supplied',...
        'GradObj','on','Display','off','HessFcn',@ObjFnHessianIPAsym,'MaxIter',15);
    problem = createOptimProblem('fmincon',...
        'lb',b_l,'ub',b_u,'x0',x0,'Aineq',A,'bineq',b,'options',options);
    [x, fTrain, ~, ~, ~ , ~ , ~] = fmincon(@ObjFnIPAsym, x0, A, b, [], [], b_l, ...
        b_u, [], options);

    % Get inputs for validation set
    [Rcut,~,~,~,~,~] = CutoffCalc2(StormVal);
    [Y, ~] = VHolCNLS(StormVal, Rcut);
    [D, DAsymC, Vm, V, StormVal] = OptDataAsym1(StormVal, Rcut);
    OptDataAsym2(DAsymC, noWvn); 
    save('Input.mat', 'Y', 'D', 'Vm', 'V', 'noWvn'); [N_Val,~] = size(D); 
    
    % Get results (validation)
    [fVal, ~] = ObjFnIPAsym(x); 
    
    % Get SSE
    SSE_train = fTrain/N_Train; SSE_val = fVal/N_Val;
end

function [SSE_Irene, SSE_Katrina, SSE_Sandy, SSE_Isabel] ...
        = ResultsTestHWind(x, noWvn, opt)

    SSE_Irene = zeros(1,5);
    SSE_Katrina = zeros(1,5);
    SSE_Sandy = zeros(1,5);
    SSE_Isabel = zeros(1,5);

    path = '/Users/changd/Dropbox (MIT)/Wind field paper (JAMC)/MATLAB_revision/';
    
    for i = 1:5
        load([path 'Training Data/TestStorms/IreneHW_' opt '.mat']);
        SSE_Irene(i) = ResultsTest(Storm, x(:,i), noWvn);
        load([path 'Training Data/TestStorms/KatrinaHW_' opt '.mat']);
        SSE_Katrina(i) = ResultsTest(Storm, x(:,i), noWvn);
        load([path 'Training Data/TestStorms/SandyHW_' opt '.mat']);
        SSE_Sandy(i) = ResultsTest(Storm, x(:,i), noWvn);
        load([path 'Training Data/TestStorms/IsabelHW_SHIPS.mat']);
        SSE_Isabel(i) = ResultsTest(Storm, x(:,i), noWvn);
    end
end

function [SSE_Katrina, SSE_Sandy] = ResultsTestWRF(x, noWvn, opt)

    SSE_Katrina = zeros(5,1);
    SSE_Sandy = zeros(5,10);
    idxSandy = [13 14 16 30 34 42 46 52 56 59];
    for i = 1:5
        %i
        load(['Training Data/KatrinaWRF_' opt '.mat']);
        SSE_Katrina(i) = ResultsTest(Storm, x(:,i), noWvn);
        for j = 1:10
            %j
            load(['Training Data/TestSandyWRF_' opt '_' num2str(idxSandy(j)) '.mat']);
            SSE_Sandy(i,j) = ResultsTest(Storm, x(:,i), noWvn);
        end
    end
end

function [SSE] = ResultsTest(Storm, xOpt, noWvn)
    [~,L]=size(xOpt);
    SSE=zeros(1,L);

    for i=1:L
        [Rcut,~,~,~,~,~] = CutoffCalc(Storm);
        
        Y = VHolCNLS(Storm,Rcut);
        
        [D,Vm,V,~,~,Storm] = OptDataAsym(Storm,Rcut, noWvn);
        
        save('Training.mat','Storm');
        save('Input.mat','Y','Rcut','D','Vm','V', 'noWvn');
        [f,df] = ObjFnIPAsym(xOpt(:,i));
        [N,~]=size(D);
        SSE(i)=f/N;
    end
end