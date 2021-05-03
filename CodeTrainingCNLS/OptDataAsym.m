 function [D,Vm,V,A,b,Storm] = OptDataAsym(Storm, Rcut, noWvn)
    
    %% Input:
    % Storm struct from previous "machine learning" scheme
    % k: Degree of polynomial
    
    %% Output:
    % x: Input/domain
    % Vm: Values for empirical max velocities (also 'input')
    % V: Output
    % A and b: Input for inequality constraints
    % Storm: Revised struct with sampled data to be used
    
    T = length(Storm);
    Vsts = zeros(1,T);
    Vshs = zeros(1,T);
    dirDiff = zeros(1,T);
    %% Produce sampled radius/velocity data
    for i=1:T      
        Vsts(i) = Storm(i).Vtrans;
        Vshs(i) = Storm(i).Vshr;
        if Storm(i).Vshrdir > 180
            Storm(i).Vshrdir = Storm(i).Vshrdir - 360;
        elseif Storm(i).Vshrdir < -180
            Storm(i).Vshrdir = Storm(i).Vshrdir + 360;
        end
        dirDiff(i) = Storm(i).Vtransdir-Storm(i).Vshrdir;
        if dirDiff(i) > 180
            dirDiff(i) = dirDiff(i) - 360;
        elseif dirDiff(i) < -180
            dirDiff(i) = dirDiff(i) + 360;
        end
        dirDiff(i) = dirDiff(i)*(pi/180);
    end    
    
    %% Initialize output
    V=[];
    VmA=[];
    Ddir=[];
    DdirSh=[];
    DAsymSt=[];
    DAsymCSt=[];
    DAsymSh=[];
    DAsymCSh=[];
    DAsymDirDiff = [];
    DAsymCDirDiff = [];

    for i=1:T
%         Vm=max(Storm(i).Vavg);
        Vm=Storm(i).Vm;
        idx=find(Storm(i).AsymR<=Rcut(i));
        V=[V Storm(i).AsymV(idx)];
        dir=Storm(i).AsymDir(idx);
        dirSh=Storm(i).AsymDirShr(idx);
        N=length(Storm(i).AsymV(idx)); 

        Ddir=[Ddir dir];
        DdirSh = [DdirSh dirSh];
        DAsymSt=[DAsymSt repmat(Vsts(i).^(0:1)',1,N)];
        DAsymCSt=[DAsymCSt Vsts(i).^(0:1)'];
        DAsymSh=[DAsymSh repmat(Vshs(i).^(0:1)',1,N)];
        DAsymCSh=[DAsymCSh Vshs(i).^(0:1)']; 
        DAsymDirDiff = [DAsymDirDiff repmat(dirDiff(i).^(0:1)',1,N)];
        DAsymCDirDiff = [DAsymCDirDiff dirDiff(i).^(0:1)'];
        VmA=[VmA Vm*ones(1,N)];
    end

    %% Produce output vectors/matrices (asymmetric component)
    D=[Ddir' DdirSh' DAsymSt' DAsymSh' DAsymDirDiff'];
    
    AASt=[-DAsymCSt' zeros(1*T,2)];
    APhiSt=[zeros(1*T,2) DAsymCSt'; 
        zeros(1*T,2) -DAsymCSt'];
    bASt=zeros(1*T,1)';
    bPhi1St=260*(pi/180)*ones(T,1)';
    bPhi2St=100*(pi/180)*ones(T,1)';
    
    AASh=[-DAsymCSh' zeros(1*T,2)];
    APhiSh=[zeros(1*T,2) DAsymCDirDiff'; 
        zeros(1*T,2) -DAsymCDirDiff'];
    bASh=zeros(1*T,1)';
    bPhi1Sh=180*(pi/180)*ones(T,1)';
    bPhi2Sh=180*(pi/180)*ones(T,1)';    
    ASh = [AASh; APhiSh];
    bSh = [bASh'; bPhi1Sh'; bPhi2Sh'];
    
    ASt = [AASt zeros(1*T,4*noWvn);
        APhiSt zeros(2*T,4*noWvn)];
    if noWvn == 1
        ASh = [zeros(3*T,4) ASh];
    elseif noWvn == 2
        ASh = [zeros(3*T,4) ASh zeros(3*T,4);
            zeros(3*T,8) ASh];
    elseif noWvn == 3
        ASh = [zeros(3*T,4) ASh zeros(3*T,8);
            zeros(3*T,8) ASh zeros(3*T,4);
            zeros(3*T,12) ASh];
    else
        disp('Error')
    end
    
    A=[ASt; ASh];
    b=[bASt'; bPhi1St'; bPhi2St'; repmat(bSh,[noWvn,1])];    
    V=V';
    Vm=VmA';