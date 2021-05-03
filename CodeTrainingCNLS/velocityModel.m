function Y = velocityModel(beta, X)
    betaTr = beta(1:4); betaTr_A = betaTr(1:2); betaTr_Phi = betaTr(3:4);
%     betaSh1 = beta(5:8); betaSh1_A = betaSh1(1:2); betaSh1_Phi = betaSh1(3:4);
%     betaSh2 = beta(9:12); betaSh2_A = betaSh2(1:2); betaSh2_Phi = betaSh2(3:4);
%     betaSh3 = beta(13:16); betaSh3_A = betaSh3(1:2); betaSh3_Phi = betaSh3(3:4);
    
    Y_Input = X(:,1); Vm = X(:,2); dirTrans = X(:,3); Vtr = X(:,4); 
%     dirShr = X(:,5); Vsh = X(:,6);
    sizeInput = size(Y_Input);
    
    A_Trans = betaTr_A(1)*ones(sizeInput) + betaTr_A(2)*Vtr; 
    Phi_Trans = betaTr_Phi(1)*ones(sizeInput) + betaTr_Phi(2)*Vtr;
%     A_Shr1 = betaSh1_A(1)*ones(sizeInput) + betaSh1_A(2)*Vsh; 
%     Phi_Shr1 = betaSh1_Phi(1)*ones(sizeInput) + betaSh1_Phi(2)*Vsh;
%     A_Shr2 = betaSh2_A(1)*ones(sizeInput) + betaSh2_A(2)*Vsh; 
%     Phi_Shr2 = betaSh2_Phi(1)*ones(sizeInput) + betaSh2_Phi(2)*Vsh;
%     A_Shr3 = betaSh3_A(1)*ones(sizeInput) + betaSh3_A(2)*Vsh; 
%     Phi_Shr3 = betaSh3_Phi(1)*ones(sizeInput) + betaSh3_Phi(2)*Vsh;
    
    velTr = A_Trans .* cos(dirTrans - Phi_Trans);
%     velSh1 = A_Shr1 .* cos(1 * (dirShr - Phi_Shr1));
%     velSh2 = A_Shr2 .* cos(2 * (dirShr - Phi_Shr2));
%     velSh3 = A_Shr3 .* cos(3 * (dirShr - Phi_Shr3));
    
    Y = Y_Input .* (Vm + velTr);
%     Y = Y_Input .* (Vm + velTr + velSh1 + velSh2 + velSh3);
end