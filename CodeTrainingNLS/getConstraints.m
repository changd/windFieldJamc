function [A, b] = getConstraints(Storm)
    load('Input.mat'); N = length(Storm);
    
    Vtr = zeros(N,1); Vsh = zeros(N,1);
    for i = 1:N
        Vtr(i) = Storm(i).Vtrans; Vsh(i) = Storm(i).Vshr;
    end
    
    D_St = [ones(N,1) Vtr];
    D_Sh = [ones(N,1) Vsh];
    
    A_A_Tr = [-D_St zeros(N,2)];
    A_Phi_Tr = [zeros(N,2) D_St; zeros(N,2) -D_St];
    A_Tr = [A_A_Tr zeros(N,4*noWvn); A_Phi_Tr zeros(2*N,4*noWvn)];
    
    A_A_Sh = [-D_Sh zeros(N,2)];
    A_Phi_Sh = [zeros(N,2) D_Sh; zeros(N,2) -D_Sh];
    A_Sh = [A_A_Sh; A_Phi_Sh]; 
    
    if noWvn == 1
        A_Sh = [zeros(3*N,4) A_Sh];
    elseif noWvn == 2
        A_Sh = [zeros(3*N,4) A_Sh zeros(3*N,4);
            zeros(3*N,8) A_Sh];
    elseif noWvn == 3
        A_Sh = [zeros(3*N,4) A_Sh zeros(3*N,8);
            zeros(3*N,8) A_Sh zeros(3*N,4);
            zeros(3*N,12) A_Sh];
    else
        disp('Error')
    end   
    A = [A_Tr; A_Sh];
    
    b_A_Tr=zeros(N,1); 
%     b_Phi_Tr1=260*(pi/180)*ones(N,1); b_Phi_Tr2=100*(pi/180)*ones(N,1);  
    b_Phi_Tr1=180*(pi/180)*ones(N,1); b_Phi_Tr2=180*(pi/180)*ones(N,1);  
    b_Tr = [b_A_Tr; b_Phi_Tr1; b_Phi_Tr2];
    
    b_A_Sh = zeros(N,1);
    b_Phi_Sh1 = 180*(pi/180)*ones(N,1); b_Phi_Sh2 = 180*(pi/180)*ones(N,1);   
    b_Sh = [b_A_Sh; b_Phi_Sh1; b_Phi_Sh2];
    
    b = [b_Tr; repmat(b_Sh, [noWvn, 1])];
end