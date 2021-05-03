function dr2_ShProd = CalcHessianSh(r, Y, D_A_St, D_Phi_St, A_St, Phi_St, dir, dirShr)
    dr2APhi_ShProd=zeros(2,2,N);
    dr2Phi2_ShProd=zeros(2,2,N);
    for i=1:N
        dr2APhi_ShProd(:,:,i)=-r(i)*D_A_Sh(i,:)'*D_Phi_Sh(i,:)*sin(dirShr(i)-Phi_Sh(i))*Y(i);
        dr2Phi2_ShProd(:,:,i)=r(i)*D_Phi_Sh(i,:)'*D_Phi_Sh(i,:)*A_Sh(i)*cos(dirShr(i)-Phi_Sh(i))*Y(i);
    end
