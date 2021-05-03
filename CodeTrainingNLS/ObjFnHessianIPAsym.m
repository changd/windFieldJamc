function df2 = ObjFnHessianIPAsym(x,lambda)
    %% Input:
    % c: Current unknown coefficient values
    % lambda: Needed for calculating second-derivatives for Lagrangian
    % components. Currently not used -- second derivatives are zero.
    
    %% Output:
    % df2: Hessian matrix

    load('Input.mat');
    
    %% Read and process input/domain
    dir=D(:,1);
    dirShr = D(:,2);
    D_A_St=D(:,3:4);
    D_Phi_St=D(:,3:4);
    D_A_Sh=D(:,5:6);
    D_Phi_Sh=D(:,7:8);
    N=length(dir);

    xSt = x(1:4);
    xSh = x(5:length(x));
    
    xA_St=xSt(1:2);
    xPhi_St=xSt(3:4); 
    
    xA_Sh = zeros(noWvn,2);
    xPhi_Sh = zeros(noWvn,2);   
    for i = 1:noWvn
        xA_Sh(i,:) = xSh(noWvn*4-3:noWvn*4-2);
        xPhi_Sh(i,:) = xSh(noWvn*4-1:noWvn*4);
    end 
    
    %% Calculate components of asymmetry
    A_St=(xA_St'*D_A_St')';
    Phi_St=(xPhi_St'*D_Phi_St')';
    A_Sh=(xA_Sh*D_A_Sh')';
    Phi_Sh=(xPhi_Sh*D_Phi_Sh')';
    VmA=Vm+A_St.*cos(dir-Phi_St);
    for i = 1:noWvn
        VmA = VmA+A_Sh(:,i).*cos(dirShr-Phi_Sh(:,i));
    end

    %% Calculate residuals
    Vpred=VmA.*Y';
    r=V-Vpred;

    %% Calculate gradient 
    J_A_St=zeros(N,2);
    J_Phi_St=zeros(N,2);
    J_A_Sh=zeros(N,2);
    J_Phi_Sh=zeros(N,2);
    J = [];
    for i=1:N
        J_A_St(i,:)=-D_A_St(i,:)*cos(dir(i)-Phi_St(i))*Y(i);
        J_Phi_St(i,:)=-A_St(i)*D_Phi_St(i,:)*sin(dir(i)-Phi_St(i))*Y(i);
    end
    J = [J_A_St J_Phi_St];
    
    for i = 1:noWvn
        for j = 1:N
            J_A_Sh(j,:)=-D_A_Sh(j,:)*cos(i*(dirShr(j)-Phi_Sh(j)))*Y(j);
            J_Phi_Sh(j,:)=-A_Sh(j)*D_Phi_Sh(j,:)*i*sin(i*(dirShr(j)-Phi_Sh(j)))*Y(j);
        end
        J = [J J_A_Sh J_Phi_Sh];
    end

    %% Calculate Hessian of asymmetry
    dr2_StProd = CalcHessianSmallMatrix(r, Y, D_A_St, D_Phi_St, A_St, Phi_St, dir, 1);  
    dr2_ShProd1 = CalcHessianSmallMatrix(r, Y, D_A_Sh, D_Phi_Sh, A_Sh, Phi_Sh, dirShr, 1);  
    dr2_ShProd2 = CalcHessianSmallMatrix(r, Y, D_A_Sh, D_Phi_Sh, A_Sh, Phi_Sh, dirShr, 2);  
    dr2_ShProd3 = CalcHessianSmallMatrix(r, Y, D_A_Sh, D_Phi_Sh, A_Sh, Phi_Sh, dirShr, 3);  

    if noWvn == 1
        dr2Prod = [dr2_StProd zeros(4,4);
              zeros(4,4) dr2_ShProd1];
    elseif noWvn == 2
        dr2Prod = [dr2_StProd zeros(4,8);
              zeros(4,4) dr2_ShProd1 zeros(4,4);
              zeros(4,8) dr2_ShProd2];          
    elseif noWvn == 3
        dr2Prod = [dr2_StProd zeros(4,12);
              zeros(4,4) dr2_ShProd1 zeros(4,8);
              zeros(4,8) dr2_ShProd2 zeros(4,4);
              zeros(4,12) dr2_ShProd3];    
    end
    df2=J'*J+dr2Prod;
