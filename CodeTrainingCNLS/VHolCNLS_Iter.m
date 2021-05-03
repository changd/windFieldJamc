function [Y, R, idx] = VHolCNLS_Iter(R, Vm, Rmax, B, Rn, Xn, Rcut)
    
    Xcalc = @(Rm, Rn, Xn, R) 0.5+(R-Rm)*(Xn-0.5)/(Rn-Rm);
    V_HRev=@(Vm,Rm,B,R,N,X) ((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^X;
    
    idx=find(R < Rcut); R = R(idx);

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