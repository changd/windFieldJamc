function [x, ci, SSE_train] = getResultsStorm(Storm, x0, b_l, b_u, noWvn, Rbounds)
    
    % Filter by latitude
    lat = zeros(length(Storm), 1);
    for i = 1:length(Storm)
        lat(i) = Storm(i).cent(1);
    end
    idx = find(lat >= 16 & lat <= 28);
    Storm = Storm(idx);
    
    % Get inputs for optimization problem (training set)    
    [X, V, R] = getInputs(Storm, Rbounds); 
    V_Hol = updateStormAsymTrans(Storm, Rbounds); N_Train = length(V);
    V = V-V_Hol;    
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    [A, b] = getConstraints(Storm); lambda=0;

    % Solve optimization problem
    velocityModel(x0); velocityModelHessian(x0,lambda);
    options = optimoptions('fmincon','Algorithm','interior-point','Hessian','user-supplied',...
        'GradObj','on','Display','off','HessFcn',@velocityModelHessian,'MaxIter',100);
    problem = createOptimProblem('fmincon',...
        'lb',b_l,'ub',b_u,'x0',x0,'Aineq',A,'bineq',b,'options',options);
    [x, fTrain, ~, ~, ~ , ~ , ~] = fmincon(@velocityModel, x0, A, b, [], [], b_l, ...
        b_u, [], options);
    
    % Get confidence interval
    [X, V, R] = getInputs(Storm, Rbounds); 
    V_Hol = updateStormAsymTrans(Storm, Rbounds); N_Train = length(V);
    V = V-V_Hol;   
    save('Input.mat', 'X', 'V', 'R', 'noWvn');
    velocityModel(x);
    load('residuals.mat');
    
    if sum(abs(x(5:end))) > 0
        ci = nlparci(x, r, 'jacobian', J);
    else
        ci = nlparci(x(1:4), r, 'jacobian', J(:,1:4));
    end
    
    % Get SSE
    SSE_train = fTrain/N_Train;
end

function [X, V, R] = getInputs(Storm, Rbounds)
    R = [];
    Y = []; Vm = []; dirTrans = []; Vtr = []; dirShr = []; Vsh = []; V = [];
    [Rcut, ~, ~, ~, ~, ~] = CutoffCalc2(Storm);
    for i = 1:length(Storm)
        [YIter, ~, idx] = VHolCNLS_Iter(Storm(i).AsymR, Storm(i).Vm, ...
            Storm(i).Rmax, Storm(i).B, Storm(i).Rn, Storm(i).Xn, Rcut(i), ...
            Rbounds);

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

function [Y, R, idx] = VHolCNLS_Iter(R, Vm, Rmax, B, Rn, Xn, Rcut, Rbounds)
    
    Xcalc = @(Rm, Rn, Xn, R) 0.5+(R-Rm)*(Xn-0.5)/(Rn-Rm);
    V_HRev=@(Vm,Rm,B,R,N,X) ((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^X;
    
    idx=find(R < Rcut & R/Rmax >= Rbounds(1) & R/Rmax <= Rbounds(2)); 
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

function [f, df] = velocityModel(beta)

    load('Input.mat');
    V_data = V;
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

function dir = getDirNormal(Storm, index)
    dir = Storm(index).AsymDir*(360/(2*pi));
    dir = dir + Storm(index).Vtransdir*ones(size(Storm(index).AsymDir));
    idx = find(dir < 0);
    dir(idx) = dir(idx) + 360*ones(size(dir(idx)));
end