function [MSE, SSE, residuals] = updateStormAsymTrans(Storm)
    
    residuals = []; N = 0;
    [Rcut, ~, ~, ~, ~, ~] = CutoffCalc2(Storm);
    for i = 1:length(Storm)
        % Get lambdas with respect to north + Holland velocities 
        dir = getDirNormal(Storm, i);
        [Y, ~, idx] = VHolCNLS_Iter(Storm(i).AsymR, Storm(i).Vm, Storm(i).Rmax, ...
            Storm(i).B, Storm(i).Rn, Storm(i).Xn, Rcut(i));
        dir = dir(idx);
        
        % Break down MF to u- and v-components
        [u, v] = getMF_Components(Y, dir);  
        
        % Break down translation vector to u- and v-components
        [uTrans, vTrans] = getTranslationComponents(Storm(i).Vtrans, ...
            Storm(i).Vtransdir);
    
        % Vector addition
        uAsym = u + uTrans*ones(size(u));
        vAsym = v + vTrans*ones(size(u));
        
        % Get velocity magnitude
        velAsym = sqrt(uAsym.^2 + vAsym.^2);
        
        % Update storm structure
        Storm(i).velAsymTrans = velAsym;
        
        residuals_Iter = (Storm(i).velAsymTrans - Storm(i).AsymV(idx)').^2;
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
    V_HRev=@(Vm,Rm,B,R,N,X) ((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^X;
    
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