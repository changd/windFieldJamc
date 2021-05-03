function [Rcut,Vcut,R15,V15,Rfrac,Vfrac] = CutoffCalc2(Storm)
    Xcalc = @(Rm, Rn, Xn, R) 0.5+(R-Rm)*(Xn-0.5)/(Rn-Rm);
    V_HRev=@(Vm,Rm,B,R,N,X) Vm*((Rm./R).^B.*exp(ones(1,N)-(Rm./R).^B)).^X;
    
    T = length(Storm);
    R15=zeros(1,T); V15=zeros(1,T); Rfrac=zeros(1,T); Vfrac=zeros(1,T);

    for i=1:T
        R=0:1:350;
        idx=find(R>Storm(i).Rmax);
        R=R(idx);
        Rmax = Storm(i).Rmax; Rn = Storm(i).Rn; Xn = Storm(i).Xn;
        Vm = Storm(i).Vm; B = Storm(i).B;
        X = zeros(size(R));
        for j = 1:length(R)
            if R(j) > Rmax
                X(j) = Xcalc(Rmax, Rn, Xn, R(j));
            else
                X(j) = 0.5;
            end
        end
        VHol = V_HRev(Vm, Rmax, B, R, length(R), X);
        [~,idx1]=min(abs(VHol-15*ones(size(VHol))));
        [~,idx2]=min(abs(VHol-(3/4)*Vm*ones(size(VHol))));
        if isempty(idx)==0
            R15(i)=R(idx1); V15(i)=VHol(idx1);
            Rfrac(i)=R(idx2); Vfrac(i)=VHol(idx2);
            Vcut(i)=min(V15(i),Vfrac(i));
            X = Xcalc(Rmax, Rn, Xn, 350);
            Vcut(i)=max(V_HRev(Vm,Rmax,B,350,1,X),Vcut(i));
            Rcut(i)=max(R15(i),Rfrac(i));
            Rcut(i)=min(Rcut(i),300);   
        else
            Rcut(i)=300;
            X = Xcalc(Rmax, Rn, Xn, 300);
            Vcut(i)=V_HRev(Vm,Rmax,B,300,1,X);
            R15(i)=Rcut(i); Rfrac(i)=Rcut(i);
            V15(i)=Vcut(i); Vfrac(i)=Vcut(i);
        end     
    end
        