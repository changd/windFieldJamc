function [Y, R] = VHolCNLS(Storm,Rcut)
    Xcalc = @(Rm, Rn, Xn, R) 0.5+(R-Rm)*(Xn-0.5)/(Rn-Rm);
    V_H=@(Vm,Rm,B,R,N) ((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^0.5;
    V_HRev=@(Vm,Rm,B,R,N,X) ((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^X;
    
    Y=[]; R = []; T = length(Storm);
    lengthYs = zeros(1, length(Storm)); lengthRs = zeros(1, length(Storm));
    for i=1:T
        idx=find(Storm(i).AsymR<Rcut(i));
        Vm = Storm(i).Vm; Rmax = Storm(i).Rmax; B = Storm(i).B;
        AsymR = Storm(i).AsymR(idx);
        Rn = Storm(i).Rn; Xn = Storm(i).Xn; 
        
        X = zeros(1, length(AsymR));
        for j = 1:length(AsymR)
            if AsymR(j) > Rmax;
                X(j) = Xcalc(Rmax, Rn, Xn, AsymR(j));
            else
                X(j) = 0.5;
            end
        end
%         VHol = V_H(Vm, Rmax, B, AsymR, length(AsymR));
        VHol=V_HRev(Vm, Rmax, B, AsymR, length(AsymR), X);
        Y=[Y VHol];
        R = [R AsymR/Storm(i).Rmax];
        
%         figure
%         plot(Storm(i).AsymR, Storm(i).AsymV, '.')
%         hold on
%         plot(AsymR, VHol, '.', 'MarkerSize', 16)
    end